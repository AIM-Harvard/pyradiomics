# This is a simple script to check if your parameter file is valid. Run it with the location of your parameter
# file as a command line argument (i.e. 'python testParams.py PATH/TO/PARAMFILE'). If successful, a message displaying
# custom parameters specified will be printed. If validation fails, an error message specifying cause of validation
# error will be printed.

import argparse

import pykwalify.core

from radiomics import getParameterValidationFiles


def main(paramsFile, is_model=False):
  if is_model:
    validate_model_file(paramsFile)
  else:
    validate_customization(paramsFile)


def validate_model_file(model_file):
  schema_data, schemaFuncs = getParameterValidationFiles(is_model_validation=True)
  c = pykwalify.core.Core(source_file=model_file, schema_data=schema_data, extensions=[schemaFuncs])

  try:
    params = c.validate()
    print('Model validation successfull!\n\n'
          '###Model Type###\n%s\n'
          % (params['model']['name']))
  except pykwalify.core.SchemaError as e:
    print('Parameter validation failed!\n%s' % e.msg)


def validate_customization(parameter_file):
  schema_data, schemaFuncs = getParameterValidationFiles()
  c = pykwalify.core.Core(source_file=parameter_file, schema_data=schema_data, extensions=[schemaFuncs])

  try:
    params = c.validate()
    print('Parameter validation successfull!\n\n'
          '###Enabled Features###\n%s\n'
          '###Enabled Image Types###\n%s\n'
          '###Settings###\n%s' % (params['featureClass'], params['imageType'], params['setting'])
          )
  except pykwalify.core.SchemaError as e:
    print('Parameter validation failed!\n%s' % e.msg)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('parameter_file', help='File representing the yaml or json structured configuration file to be '
                                             'tested')
  parser.add_argument('--model', '-m', action='store_true',
                      help='If this argument is specified, the configuration file is treated as a PyRadiomics Model, '
                           'otherwise, it is treated as an extraction parameter file')

  args = parser.parse_args()
  main(args.parameter_file, args.model)
