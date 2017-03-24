This folder contains the test cases and baseline used to ensure consistent results from features implemented in 
PyRadioimcs.

`baseline` contains calculated features in full-python mode for test cases brain1, brain2, breast1, lung1 and lung2 and
was generated using `addClassToBaseline.py`, available in the `bin` folder of this repository.

Settings used were hardcoded in the script to generate the baseline and apply no filters or resampling.
Applied settings are contained within the general_info section included in the baseline. This
section must be present, as the test code extracts the settings to apply from the baseline, enabling
a baseline to have different settings if necessary.
