# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_exampleSettings.py

import os

from nose_parameterized import parameterized
import pykwalify.core

from radiomics import getParameterValidationFiles


def exampleSettings_name_func(testcase_func, param_num, param):
  return '%s_%s' % (testcase_func.__name__, os.path.splitext(os.path.basename(param.args[0]))[0])


class TestExampleSettings:
  def __init__(self):
    self.schemaFile, self.schemaFuncs = getParameterValidationFiles()

  def generateScenarios():
    dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'examples', 'exampleSettings')
    if os.path.isdir(dataDir):
      settingsFiles = [fname for fname in os.listdir(dataDir) if fname.endswith('.yaml') or fname.endswith('.yml')]

      for fname in settingsFiles:
        yield os.path.join(dataDir, fname)

  @parameterized.expand(generateScenarios(), testcase_func_name=exampleSettings_name_func)
  def test_scenarios(self, settingsFile):

    assert os.path.isfile(self.schemaFile)
    assert os.path.isfile(self.schemaFuncs)
    c = pykwalify.core.Core(source_file=settingsFile, schema_files=[self.schemaFile], extensions=[self.schemaFuncs])
    c.validate()
