import argparse
import logging

import six

import radiomics
from radiomics import featureextractor
from testUtils import RadiomicsTestUtils


def main(argv=None):
  logger = logging.getLogger('radiomics.addTest')

  logger.setLevel(logging.INFO)
  radiomics.setVerbosity(logging.INFO)

  handler = logging.FileHandler(filename='testLog.txt', mode='w')
  formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
  handler.setFormatter(formatter)

  radiomics.logger.addHandler(handler)

  parser = argparse.ArgumentParser()
  parser.add_argument('TestName', type=str, help='Name for the new test, must not be already present in the baseline')
  parser.add_argument('TestCase', type=str, choices=list(radiomics.testCases), help='Test image and segmentation to '
                                                                                    'use in the new test')
  parser.add_argument('Configuration', metavar='FILE', default=None,
                      help='Parameter file containing the settings to be used in extraction')
  parser.add_argument('--force', '-f', action='store_true', help='When the test is already known for the class, '
                                                                 'overwrite it in stead of just skipping that class')

  args = parser.parse_args(argv)

  logger.info('Input accepted, starting addTest...')

  try:

    testutils = RadiomicsTestUtils()

    try:
      assert args.TestCase in radiomics.testCases
    except AssertionError:
      logger.error('Input not valid, cancelling addTest!')
      exit(1)

    logger.debug('Initializing extractor')
    extractor = featureextractor.RadiomicsFeaturesExtractor(args.Configuration)

    logger.debug('Starting extraction')
    featurevector = extractor.execute(*radiomics.getTestCase(args.TestCase))

    configuration = {}
    baselines = {}

    for k, v in six.iteritems(featurevector):
      if 'diagnostics' in k:
        configuration[k] = v
      else:
        image_filter, feature_class, feature_name = k.split('_')
        if feature_class not in baselines:
          baselines[feature_class] = {}
        baselines[feature_class][k] = v

    configuration['diagnostics_Configuration_TestCase'] = args.TestCase

    testutils.addTest(args.TestName, configuration, baselines)

    logger.info('addTest Done')

  except Exception:
    logger.error('Error running addTest!', exc_info=True)


if __name__ == '__main__':
  main()
