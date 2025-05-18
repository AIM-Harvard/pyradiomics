from __future__ import annotations

import logging

from radiomics import getFeatureClasses

featureClasses = getFeatureClasses()


def pytest_generate_tests(metafunc):
    metafunc.parametrize(
        ["featureClassName", "featureName"], metafunc.cls.generate_scenarios()
    )


class TestDocStrings:

    @staticmethod
    def generate_scenarios():
        for featureClassName, featureClass in featureClasses.items():
            logging.info("generate_scenarios %s", featureClassName)
            doc = featureClass.__doc__
            assert doc is not None

            featureNames = featureClass.getFeatureNames()
            for f in featureNames:
                yield (featureClassName, f)

    def test_class(self, featureClassName, featureName):
        logging.info("%s", featureName)
        features = featureClasses[featureClassName]
        doc = getattr(features, f"get{featureName}FeatureValue").__doc__
        logging.info("%s", doc)
        assert doc is not None
