# This is a simple script to check if your parameter file is valid. Run it with the location of your parameter
# file as a command line argument (i.e. 'python testParams.py PATH/TO/PARAMFILE'). If successful, a message displaying
# custom parameters specified will be printed. If validation fails, an error message specifying cause of validation
# error will be printed.

import sys

import pykwalify.core

from radiomics import getParameterValidationFiles

def main(paramsFile):
  schemaFile, schemaFuncs = getParameterValidationFiles()

  c = pykwalify.core.Core(source_file=paramsFile, schema_files=[schemaFile], extensions=[schemaFuncs])
  try:
    params = c.validate()
    print('Parameter validation successfull!\n\n'
          '###Enabled Features###\n%s\n'
          '###Enabled Image Types###\n%s\n'
          '###Settings###\n%s' % (params['featureClass'], params['imageType'], params['setting']))
  except Exception as e:
    print('Parameter validation failed!\n%s' % e.message)


if __name__ == '__main__' and len(sys.argv) > 1:
    main(sys.argv[1])
