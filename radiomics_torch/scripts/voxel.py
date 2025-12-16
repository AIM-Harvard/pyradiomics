from __future__ import annotations

import logging.config
import os
import threading
from collections import OrderedDict
from datetime import datetime

import SimpleITK as sitk

import radiomics_torch

caseLogger = logging.getLogger("radiomics.script")


class _SingletonVoxelParallelExtractionConfig:
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
                    _SingletonVoxelParallelExtractionConfig._initialized = True


def extractVoxel(case_idx, case, extractor, **kwargs):
    out_dir = kwargs.get("out_dir")
    unix_path = kwargs.get("unix_path", False)

    # Instantiate the output
    feature_vector = OrderedDict(case)

    try:
        if out_dir is None:
            out_dir = "."
        elif not os.path.isdir(out_dir):
            caseLogger.debug(f"Creating output directory at {out_dir}")
            os.makedirs(out_dir)

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
        result = extractor.execute(
            imageFilepath, maskFilepath, label, label_channel, voxelBased=True
        )

        for k in result:
            if isinstance(result[k], sitk.Image):
                target = os.path.join(out_dir, f"Case-{int(case_idx)}_{k}.nrrd")
                sitk.WriteImage(result[k], target, True)
                if unix_path and os.path.sep != "/":
                    target = target.replace(os.path.sep, "/")
                feature_vector[k] = target
            else:
                feature_vector[k] = result[k]

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
        caseLogger.error(f"Feature extraction failed! {e}", exc_info=True)
        # Log but don't raise error

    return feature_vector


def extractVoxel_parallel(args, logging_config=None, **kwargs):
    try:
        # set thread name to patient name
        threading.current_thread().name = f"case {args[0]}"  # args[0] = case_idx

        if logging_config is not None:
            _configureParallelExtraction(logging_config)

        return extractVoxel(*args, **kwargs)
    except (KeyboardInterrupt, SystemExit):
        # Catch the error here, as this represents the interrupt of the child process.
        # The main process is also interrupted, and cancellation is further handled there
        return None


def _configureParallelExtraction(logging_config, add_info_filter=True):
    """
    Initialize logging for parallel extraction. This needs to be done here, as it needs to be done for each thread that is
    created.
    """

    if _SingletonVoxelParallelExtractionConfig().parallel_extraction_configured:
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

    _SingletonVoxelParallelExtractionConfig().parallel_extraction_configured = True
    radiomics_torch.logger.debug("parallel extraction configured")
