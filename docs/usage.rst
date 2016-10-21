=====
Usage
=====

-------
Example
-------

* Run the helloRadiomics example, using sample data provided in ``pyradiomics/data``:

  * ``python bin/helloRadiomics.py``

---------------
Interactive Use
---------------

* Add pyradiomics to the environment variable PYTHONPATH:

  *  ``setenv PYTHONPATH /path/to/pyradiomics/radiomics``

* Start the python interactive session:

  * ``python``

* Import the necessary classes::

     from radiomics import firstorder, glcm, imageoperations, shape, glrlm, glszm
     import SimpleITK as sitk
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
     for (key,val) in firstOrderFeatures.featureValues.iteritems():
       print '  ',key,':',val

* See the :ref:`radomics package<radiomics-firstorder-label>` for more features that you can calculate.
