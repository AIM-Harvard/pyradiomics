.. _radiomics-contributing-label:

===========================
Contributing to pyradiomics
===========================

There are many ways to contribute to pyradiomics, with varying levels of effort.  Do try to
look through the documentation first if something is unclear, and let us know how we can
do better.

* Ask a question on the `pyradiomics email list <https://groups.google.com/forum/#!forum/pyradiomics>`_
* Submit a parameter file you used for your extraction
* Submit a feature request or bug, or add to the discussion on the `pyradiomics issue tracker <https://github.com/Radiomics/pyradiomics/issues>`_
* Submit a `Pull Request <https://github.com/Radiomics/pyradiomics/pulls>`_ to improve pyradiomics or its documentation

We encourage a range of Pull Requests, from patches that include passing tests and
documentation, all the way down to half-baked ideas that launch discussions.

----------------------------------------------
The PR Process, Circle CI, and Related Gotchas
----------------------------------------------

How to submit a PR ?
####################

If you are new to pyradiomics development and you don't have push access to the pyradiomics
repository, here are the steps:

1. `Fork and clone <https://help.github.com/articles/fork-a-repo/>`_ the repository.
2. Create a branch.
3. `Push <https://help.github.com/articles/pushing-to-a-remote/>`_ the branch to your GitHub fork.
4. Create a `Pull Request <https://github.com/Radiomics/pyradiomics/pulls>`_.

This corresponds to the *Fork & Pull Model* mentioned in the `GitHub flow <https://guides.github.com/introduction/flow/index.html>`_
guides.

If you have push access to pyradiomics repository, you could simply push your branch
into the main repository and create a `Pull Request <https://github.com/Radiomics/pyradiomics/pulls>`_. This corresponds to the
*Shared Repository Model* and will facilitate other developers to checkout your
topic without having to `configure a remote <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_.
It will also simplify the workflow when you are *co-developing* a branch.

When submitting a PR, make sure to add a ``cc: @Radiomics/developers`` comment to notify pyradiomics
developers of your awesome contributions. Based on the
comments posted by the reviewers, you may have to revisit your patches.

How to integrate a PR ?
#######################

Getting your contributions integrated is relatively straightforward, here
is the checklist:

* Your changes include an update of the documentation if necessary

  * Documentation on modules, classes and functions is contained in the respective docstrings
  * More global documentation is contained in the ``docs`` folder.
  * New modules need to be added to the auto-generated documentation. See
    `here <http://pyradiomics.readthedocs.io/en/latest/developers.html#documentation>`_ for more
    information on adding new modules to the documentation.

* Your changes are added in the `changelog <https://github.com/Radiomics/pyradiomics/tree/master/CHANGES.rst>`_ in the *Next Release*  section.
* All tests pass
* Consensus is reached. This usually means that at least one reviewer reviewed and approved your
  changes or added a ``LGTM`` comment, which is an acronym for *Looks Good to Me*.

Next, there are two scenarios:

* You do NOT have push access: A pyradiomics core developer will integrate your PR.
* You have push access: Simply click on the "Merge pull request" button.

Then, click on the "Delete branch" button that appears afterward.

Automatic testing of pull requests
##################################

Every pull request is tested automatically using CircleCI, TravisCI and AppVeyor each time you push
a commit to it. The Github UI will restrict users from merging pull requests until
the builds have returned with a successful result indicating that all tests have
passed and there were no problems detected by the linter. These tests include the following

* flake8 to check adherence to the code style. See ``.flake8`` and ``.editorconfig`` for styles,
  exceptions to the PEP8 style, etc.
* If a feature class has a function ``_calculateCMatrix()``, identifying it as a C enhanced class,
  output from the C extension is compared to the output from full python calculation. A absolute
  difference of 1e-3 is allowed to account for machine precision errors.
* All implemented features and feature classes have docstrings at the class level and feature
  definition level.
* A baseline is available for all features extracted from the 5 included test cases and
  calculated features match this baseline to within 3% (allowing for machine precision errors)

.. _radiomics-submit-parameter-file-label:

---------------------------
Submitting a parameter file
---------------------------

Different inputs into PyRadiomics require different settings. We encourage users to share their parameter file to help
others extract features using the best settings for their use case.

How to submit your parameter file?
##################################

Parameter files are stored in the repository under ``examples/exampleSettings``. If you wish to submit your parameters
to the community, you can add your file here via a pull request (see above for details on making PRs).
To help you along, here is a small checklist:

* The filename should at least contain the modality (e.g. "MR") for which it is intended, and optionally the body part
  (e.g. "prostate").
* Ensure the file has the correct extension (either ".yml" or ".yaml")
* Using comments in the parameter file, briefly explain your use case.

After you've opened a PR to submit your parameter file, it will be checked for validity by the automatic testing. You
don't have to specify your file anywhere, as parameter files are detected automatically in the ``exampleSettings``
folder. If you want, you can also check your file manually using ``bin/testParams.py``.
