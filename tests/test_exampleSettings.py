from __future__ import annotations

import os

import pykwalify.core

from radiomics import getParameterValidationFiles

schemaFile, schemaFuncs = getParameterValidationFiles()


def pytest_generate_tests(metafunc):
    metafunc.parametrize("settingsFile", metafunc.cls.generate_scenarios())


def exampleSettings_name_func(testcase_func, param):
    return f"{testcase_func.__name__}_{os.path.splitext(os.path.basename(param.args[0]))[0]}"


class TestExampleSettings:
    @staticmethod
    def generate_scenarios():
        dataDir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "examples",
            "exampleSettings",
        )
        if os.path.isdir(dataDir):
            settingsFiles = [
                fname
                for fname in os.listdir(dataDir)
                if fname.endswith((".yaml", ".yml"))
            ]

            for fname in settingsFiles:
                yield os.path.join(dataDir, fname)

    def test_scenarios(self, settingsFile):
        assert os.path.isfile(schemaFile)
        assert os.path.isfile(schemaFuncs)
        c = pykwalify.core.Core(
            source_file=settingsFile,
            schema_files=[schemaFile],
            extensions=[schemaFuncs],
        )
        c.validate()
