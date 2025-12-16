from __future__ import annotations

import csv
import logging.config
import os
import threading
from collections import OrderedDict
from datetime import datetime

import SimpleITK as sitk

import radiomics_torch

caseLogger = logging.getLogger("radiomics.script")


class _SingletonSegmentParallelExtractionConfig:
    _instance = None
    _initialized = False
    _lock = threading.Lock()

    def __new__(cls, *_args, **_kwargs):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, configured=False):
        if not self._initialized:
            with self._lock:
                if not self._initialized:
                    self.parallel_extraction_configured = configured
                    _SingletonSegmentParallelExtractionConfig._initialized = True


def extractSegment(case_idx, case, extractor, **kwargs):
    out_dir = kwargs.get("out_dir")

    if out_dir is None:
        return _extractFeatures(case_idx, case, extractor)

    filename = os.path.join(out_dir, f"features_{case_idx}.csv")
    if os.path.isfile(filename):
        # Output already generated, load result (prevents re-extraction in case of interrupted process)
        with open(filename) as outputFile:
            reader = csv.reader(outputFile)
            headers = next(reader)
            values = next(reader)
            feature_vector = OrderedDict(zip(headers, values))

        caseLogger.info("Patient %s already processed, reading results...", case_idx)
    else:
        # Extract the set of features. Set parallel_config flag to None, as any logging initialization is already handled.
        feature_vector = _extractFeatures(case_idx, case, extractor)

        # Store results in temporary separate files to prevent write conflicts
        # This allows for the extraction to be interrupted. Upon restarting, already processed cases are found in the
        # TEMP_DIR directory and loaded instead of re-extracted
        with open(filename, "w") as outputFile:
            writer = csv.DictWriter(
                outputFile, fieldnames=list(feature_vector.keys()), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerow(feature_vector)

    return feature_vector


def _extractFeatures(case_idx, case, extractor):
    # Instantiate the output
    feature_vector = OrderedDict(case)

    try:
        caseLogger.info("Processing case %s", case_idx)
        t = datetime.now()

        imageFilepath = case["Image"]  # Required
        maskFilepath = case["Mask"]  # Required
        label = case.get("Label", None)  # Optional
        if isinstance(label, str):
            label = int(label)
        label_channel = case.get("Label_channel", None)  # Optional
        if isinstance(label_channel, str):
            label_channel = int(label_channel)

        # Extract features
        feature_vector.update(
            extractor.execute(imageFilepath, maskFilepath, label, label_channel)
        )

        # Display message
        delta_t = datetime.now() - t
        caseLogger.info("Case %s processed in %s", case_idx, delta_t)

    except (
        KeyboardInterrupt,
        SystemExit,
    ):  # Cancel extraction by forwarding this 'error'
        raise
    except SystemError as e:
        # Occurs when Keyboard Interrupt is caught while the thread is processing a SimpleITK call
        raise KeyboardInterrupt() from e
    except Exception as e:
        caseLogger.error(f"Feature extraction failed! : {e}", exc_info=True)
        # log but do not raise exception

    return feature_vector


def extractSegment_parallel(args, logging_config=None, **kwargs):
    try:
        # set thread name to patient name
        threading.current_thread().name = f"case {args[0]}"  # args[0] = case_idx

        if logging_config is not None:
            _configureParallelExtraction(logging_config)

        return extractSegment(*args, **kwargs)
    except (KeyboardInterrupt, SystemExit):
        # Catch the error here, as this represents the interrupt of the child process.
        # The main process is also interrupted, and cancellation is further handled there
        return None


def _configureParallelExtraction(logging_config, add_info_filter=True):
    """
    Initialize logging for parallel extraction. This needs to be done here, as it needs to be done for each thread that is
    created.
    """
    if _SingletonSegmentParallelExtractionConfig().parallel_extraction_configured:
        return

    # Configure logging
    ###################

    logging.config.dictConfig(logging_config)

    if add_info_filter:
        # Define filter that allows messages from specified filter and level INFO and up, and level WARNING and up from
        # other loggers.
        class info_filter(logging.Filter):
            def __init__(self, name):
                super().__init__(name)
                self.level = logging.WARNING

            def filter(self, record):
                if record.levelno >= self.level:
                    return True
                return bool(record.name == self.name and record.levelno >= logging.INFO)

        # Adding the filter to the first handler of the radiomics logger limits the info messages on the output to just
        # those from radiomics.script, but warnings and errors from the entire library are also printed to the output.
        # This does not affect the amount of logging stored in the log file.
        outputhandler = radiomics_torch.logger.handlers[0]  # Handler printing to the output
        outputhandler.addFilter(info_filter("radiomics.script"))

    # Ensure the entire extraction for each cases is handled on 1 thread
    ####################################################################

    sitk.ProcessObject_SetGlobalDefaultNumberOfThreads(1)

    _SingletonSegmentParallelExtractionConfig().parallel_extraction_configured = True
    radiomics_torch.logger.debug("parallel extraction configured")
