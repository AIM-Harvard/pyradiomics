from __future__ import annotations

import logging
import os

import numpy as np
from testUtils import RadiomicsTestUtils

from radiomics import getFeatureClasses, testCases

testUtils = RadiomicsTestUtils()

featureClasses = getFeatureClasses()


def pytest_generate_tests(metafunc):
    metafunc.parametrize(
        ["testCase", "featureClassName"], metafunc.cls.generate_scenarios()
    )


class TestMatrices:

    @staticmethod
    def generate_scenarios():
        for testCase in testCases:
            if testCase.startswith("test"):
                continue
            for className, featureClass in featureClasses.items():
                assert featureClass is not None
                if "_calculateMatrix" in dir(featureClass):
                    logging.debug("generate_scenarios: featureClass = %s", className)
                    yield testCase, className

    def test_scenario(self, testCase, featureClassName):
        logging.debug(
            "test_scenario: testCase = %s, featureClassName = %s",
            testCase,
            featureClassName,
        )

        baselineFile = os.path.join(
            testUtils.getDataDir(),
            "baseline",
            f"{testCase}_{featureClassName}.npy",
        )
        assert os.path.isfile(baselineFile)

        baselineMatrix = np.load(baselineFile)

        testUtils.setFeatureClassAndTestCase(featureClassName, testCase)

        testImage = testUtils.getImage("original")
        testMask = testUtils.getMask("original")

        featureClass = featureClasses[featureClassName](
            testImage, testMask, **testUtils.getSettings()
        )
        featureClass._initCalculation()

        cMat = getattr(featureClass, f"P_{featureClassName}")
        assert cMat is not None

        # Check if the calculated arrays match
        assert np.max(np.abs(baselineMatrix - cMat)) < 1e-3
