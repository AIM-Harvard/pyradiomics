from __future__ import annotations

import logging
import os

from testUtils import RadiomicsTestUtils

from radiomics import getFeatureClasses

testUtils = RadiomicsTestUtils()
tests = sorted(testUtils.getTests())

featureClasses = getFeatureClasses()


def pytest_generate_tests(metafunc):
    metafunc.parametrize(["testCase", "featureName"], metafunc.cls.generate_scenarios())


class TestFeatures:
    @staticmethod
    def generate_scenarios():
        for test in tests:
            for featureClassName in featureClasses:
                # Get all feature names for which there is a baseline with current test case
                # Raises an assertion error when the class is not yet present in the baseline
                # Returns None if no baseline is present for this specific test case
                # Returns a list of feature names for which baseline values are present for this test
                baselineFeatureNames = testUtils.getFeatureNames(featureClassName, test)

                if baselineFeatureNames is None:
                    continue
                assert len(baselineFeatureNames) > 0

                uniqueFeatures = {f.split("_")[-1] for f in baselineFeatureNames}

                # Get a list of all features for current class
                featureNames = featureClasses[featureClassName].getFeatureNames()
                # Get a list of all non-deprecated features
                activeFeatures = {
                    f for (f, deprecated) in featureNames.items() if not deprecated
                }
                # Check if all active features have a baseline (exclude deprecated features from this set)
                if len(activeFeatures - uniqueFeatures) > 0:
                    msg = f"Missing baseline for active features {activeFeatures - uniqueFeatures}"
                    raise AssertionError(msg)
                if len(uniqueFeatures - activeFeatures) > 0:
                    msg = f"Missing function(s) for baseline feature(s) {uniqueFeatures - activeFeatures}"
                    raise AssertionError(msg)

                logging.debug(
                    "generate_scenarios: featureNames = %s", baselineFeatureNames
                )
                for featureName in baselineFeatureNames:
                    yield test, featureName

    def test_scenario(self, testCase, featureName):
        print("")
        featureClass = None

        featureName = featureName.split("_")

        logging.debug(
            "test_scenario: test = %s, featureClassName = %s, featureName = %s",
            testCase,
            featureName[1],
            featureName[-1],
        )

        testOrClassChanged = testUtils.setFeatureClassAndTestCase(
            featureName[1], testCase
        )

        testImage = testUtils.getImage(featureName[0])
        testMask = testUtils.getMask(featureName[0])

        if featureClass is None or testOrClassChanged:
            msg = f"Init {featureName[1]}"
            logging.debug(msg)
            featureClass = featureClasses[featureName[1]](
                testImage, testMask, **testUtils.getSettings()
            )

        assert featureClass is not None

        featureClass.disableAllFeatures()
        featureClass.enableFeatureByName(featureName[-1])
        featureClass.execute()
        # get the result and test it
        val = featureClass.featureValues[featureName[-1]]
        testUtils.checkResult(featureName, val)


def teardown_module():
    print("")
    res = testUtils.getResults()
    print("Results:")
    print(res)
    resultsFile = os.path.join(testUtils.getDataDir(), "PyradiomicsFeatures.csv")
    testUtils.writeCSV(res, resultsFile)
    diff = testUtils.getDiffs()
    print("Differences from baseline:")
    print(diff)
    diffFile = os.path.join(
        testUtils.getDataDir(), "Baseline2PyradiomicsFeaturesDiff.csv"
    )
    testUtils.writeCSV(diff, diffFile)
    logging.info(
        "Wrote calculated features to %s, and the differences between the baseline features and the calculated ones to %s.",
        resultsFile,
        diffFile,
    )
