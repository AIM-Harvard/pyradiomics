.. _radiomics-usage-label:

=====
Usage
=====

-----------------
Instruction Video
-----------------

.. raw:: html

        <div data-video = "ZF1rTtRW3eQ"
        data-startseconds = "221"
        data-endseconds = "538"
        data-height = "315"
        data-width = "560"
        id = "youtube-player">
        </div>

        <script src="https://www.youtube.com/iframe_api"></script>
        <script type="text/javascript">
          function onYouTubeIframeAPIReady() {
            var ctrlq = document.getElementById("youtube-player");
            var player = new YT.Player('youtube-player', {
              height: ctrlq.dataset.height,
              width: ctrlq.dataset.width,
              events: {
                'onReady': function(e) {
                  e.target.cueVideoById({
                    videoId: ctrlq.dataset.video,
                    startSeconds: ctrlq.dataset.startseconds,
                    endSeconds: ctrlq.dataset.endseconds
                  });
                }
              }
            });
          }
        </script>

-------
Example
-------

* PyRadiomics example code and data is available in the `Github repository <https://github.com/Radiomics/pyradiomics>`_

* The sample sample data is provided in ``pyradiomics/data``

* Use `jupyter <http://jupyter.org/>`_ to run the helloRadiomics example, located in ``pyradiomics/bin/Notebooks``

* Jupyter can also be used to run the example notebook as shown in the instruction video

  * The example notebook can be found in ``pyradiomics/bin/Notebooks``

  * The parameter file used in the instruction video is available in ``pyradiomics/bin``

* If jupyter is not installed, run the python script alternative (``pyradiomics/bin/helloRadiomics.py``):

  * ``python helloRadiomics.py``

----------------
Command Line Use
----------------

* PyRadiomics has 2 commandline scripts, ``pyradiomics`` is for single image feature extraction and ``pyradiomicsbatch``
  is for feature extraction from a batch of images and segmentations.

* Both scripts can be run directly from a command line window, anywhere in your system.

* To extract features from a single image and segmentation run::

    pyradiomics <path/to/image> <path/to/segmentation>

* To extract features from a batch run::

    pyradiomicsbatch <path/to/input> <path/to/output>

* The input file for batch processing is a CSV file where each row represents one combination of an image and a
  segmentation and contains 5 elements: 1) patient ID, 2) sequence name (image identifier), 3) reader (segmentation
  identifier), 4) path/to/image, 5) path/to/mask.

* For more information on passing parameter files, setting up logging and controlling output format, run::

    pyradiomics -h
    pyradiomicsbatch -h


---------------
Interactive Use
---------------

* (LINUX) Add pyradiomics to the environment variable PYTHONPATH:

  *  ``setenv PYTHONPATH /path/to/pyradiomics/radiomics``

* Start the python interactive session:

  * ``python``

* Import the necessary classes::

     from radiomics import featureextractor
     import six
     import sys, os

* Set up a pyradiomics directory variable::

    dataDir = '/path/to/pyradiomics'

* You will find sample data files brain1_image.nrrd and brain1_label.nrrd in that directory.

* Store the path of your image and mask in two variables::

    imageName = os.path.join(dataDir, "data", 'brain1_image.nrrd')
    maskName = os.path.join(dataDir, "data",  'brain1_label.nrrd')

* Also store the path to the file containing the extraction settings::

    params = os.path.join(dataDir, "bin", "Params.yaml")

* Instantiate the feature extractor class with the parameter file::

    extractor = featureextractor.RadiomicsFeaturesExtractor(params)

* Calculate the features::

    result = extractor.execute(imageName, maskName)
    for key, val in six.iteritems(result):
      print("\t%s: %s" %(key, val))

* See the :ref:`feature extractor class<radiomics-featureextractor-label>` for more information on using this core class.

------------------------
PyRadiomics in 3D Slicer
------------------------

A convenient front-end interface is provided as the 'Radiomics' extension for 3D Slicer. It is available
`here <https://github.com/Radiomics/SlicerRadiomics>`_.

------------------------------
Using feature classes directly
------------------------------

* This represents an example where feature classes are used directly, circumventing checks and preprocessing done by
  the radiomics feature extractor class, and is not intended as standard use example.

* (LINUX) Add pyradiomics to the environment variable PYTHONPATH:

  *  ``setenv PYTHONPATH /path/to/pyradiomics/radiomics``

* Start the python interactive session:

  * ``python``

* Import the necessary classes::

     from radiomics import firstorder, glcm, imageoperations, shape, glrlm, glszm
     import SimpleITK as sitk
     import six
     import sys, os

* Set up a data directory variable::

    dataDir = '/path/to/pyradiomics/data'

* You will find sample data files brain1_image.nrrd and brain1_label.nrrd in that directory.

* Use SimpleITK to read a the brain image and mask::

     imageName = str(dataDir + os.path.sep + 'brain1_image.nrrd')
     maskName = str(dataDir + os.path.sep + 'brain1_label.nrrd')
     image = sitk.ReadImage(imageName)
     mask = sitk.ReadImage(maskName)

* Calculate the first order features::

     firstOrderFeatures = firstorder.RadiomicsFirstOrder(image,mask)
     firstOrderFeatures.calculateFeatures()
     for (key,val) in six.iteritems(firstOrderFeatures.featureValues):
       print("\t%s: %s" % (key, val))

* See the :ref:`radiomics-features-label` section for more features that you can calculate.
