# About

This is an experimental script to support the use of pyradiomics with DICOM data.

The script will accept as input a directory with a single DICOM image study for the input image,
and the file name pointing to a DICOM Segmentation Image (DICOM SEG) object.

The script will transparently convert the DICOM image into a representation suitable by pyradiomics
using either plastimatch or dcm2niix.

# Why?

* medical image data usually comes in DICOM, and pyradiomics users often ask for help working with DICOM data
* there are public collections of data on TCIA where segmentations are stored as DICOM SEG
* the use of DICOM representation for radiomics features
 * introduces standardized formalism for the attributes that should be stored to accompany the features
 * allows to link results of calculations with the various ontologies describing the anatomy of the regions
  analyzed, and the features itself (e.g., the SR document produced by the script will utilize IBSI nomenclature
  to describe those features implemented in pyradiomics that have correspondence in IBSI)
 * allows to reference (by unique identifiers) the DICOM image series and DICOM segmentation used for feature 
  calculation
 * enables harmonized representation of data for images, segmentations and features (i.e., same data management
  system can be used for all data types)
 * does not prevent the use of the results in software tools that are not DICOM-aware - dcmqi can be used to 
  convert DICOM segmentations and DICOM SR with the measurements into non-DICOM representation (ITK-readable 
  image formats for segmentations, and JSON for measurements); a separate tool is available to generate 
  tab-delimited representation for DICOM attributes and measurements stored in those SRs: https://github.com/QIICR/dcm2tables

# Prerequisites

* [plastimatch](http://plastimatch.org/plastimatch.html) or [dcm2niix](https://github.com/rordenlab/dcm2niix) for image volume reconstruction
* dcmqi (build from https://github.com/QIICR/dcmqi/commit/3638930723bf1a239515409c1f9ec886a9fedb41 or later) for reading DICOM SEG and converting to a representation suitable by pyradiomics, and for storing the resulting features as a DICOM Structured Report, instantiating SR TID 1500
* prior to using this script, you might want to sort your DICOM data such that individual series
are stored in separate directories. You might find this tool useful for this purpose: https://github.com/pieper/dicomsort
* if you segmentations are not stored as DICOM SEG, you can use dcmqi for generating standard representation
of those segmentations: https://github.com/QIICR/dcmqi

# Usage

```bash
$ python pyradiomics-dcm.py -h                                                                              2.3.6
usage: pyradiomics-dcm.py --input-image <dir> --input-seg <name> --output-sr <name>

Warning: This is a "pyradiomics labs" script, which means it is an experimental feature in development!
The intent of this helper script is to enable pyradiomics feature extraction directly from/to DICOM data.
The segmentation defining the region of interest must be defined as a DICOM Segmentation image.
Support for DICOM Radiotherapy Structure Sets for defining region of interest may be added in the future.

optional arguments:
  -h, --help            show this help message and exit
  --input-image-dir <folder>
                        Path to the directory with the input DICOM series. It is expected that a single
                        series is corresponding to a single scalar volume.
  --input-seg-file <file>
                        Path to the input segmentation defined as a DICOM Segmentation object.
  --output-dir <folder>
                        Path to the directory for saving the resulting DICOM file.
  --parameters <parameters>
                        Pyradiomics feature extractor positional arguments
  --temp-dir <folder>   Path to the directory to store intermediate results
  --features-dict <file>
                        Path to the dictionary mapping pyradiomics feature names to the IBSI defined
                        features.
  --volume-reconstructor <plastimatch or dcm2niix>
                        Choose the tool to be used for reconstructing image volume from the DICOM image
                        series. Allowed options are plastimatch or dcm2niix (should be installed on the
                        system). plastimatch will be used by default.
  --geometry-tolerance <number>
                        Decimal number setting geometry tolerance for the extractor. Defaults to 1e-6.
  --correct-mask        Boolean flag argument. If present, PyRadiomics will attempt to resample the mask
                        to the image geometry if the mask check fails.
```

# Sample invocation

```bash
$ python pyradiomics-dcm.py --input-image-dir CT --input-seg SEG/1.dcm \
   --output-dir OutputSR --temp-dir TempDir --parameters Pyradiomics_Params.yaml
dcmqi repository URL: https://github.com/QIICR/dcmqi.git revision: 3638930 tag: latest-4-g3638930
Row direction: 1 0 0
Col direction: 0 1 0
Z direction: 0 0 1
Total frames: 177
Total frames with unique IPP: 177
Total overlapping frames: 0
Origin: [-227.475, -194.775, -1223]
dcmqi repository URL: https://github.com/QIICR/dcmqi.git revision: 3638930 tag: latest-4-g3638930
Total measurement groups: 1
Adding to compositeContext: 1.dcm
Composite Context initialized
SR saved!

$ dsrdump OutputSR/1.2.276.0.7230010.3.1.4.0.60427.1539113881.935517.dcm
Enhanced SR Document

Patient             : interobs05 (#interobs05)
ENH: include pyradiomics identification and version
Study               : interobs05_20170910_CT
Series              : GTV segmentation - Reader AB - pyradiomics features (#1)
Manufacturer        : QIICR (https://github.com/QIICR/dcmqi.git, #0)
Completion Flag     : PARTIAL
Verification Flag   : UNVERIFIED
Content Date/Time   : 2018-10-09 15:38:01

<CONTAINER:(,,"Imaging Measurement Report")=SEPARATE>
  <has concept mod CODE:(,,"Language of Content Item and Descendants")=(eng,RFC5646,"English")>
  <has obs context CODE:(,,"Observer Type")=(121007,DCM,"Device")>
  <has obs context UIDREF:(,,"Device Observer UID")="1.3.6.1.4.1.43046.3.1.4.0.60427.1539113880.935515">
  <has obs context TEXT:(,,"Device Observer Name")="pyradiomics">
  <has obs context TEXT:(,,"Device Observer Model Name")="2.1.0.post10.dev0+g51bc87f">
  <has concept mod CODE:(,,"Procedure reported")=(P0-0099A,SRT,"Imaging procedure")>
  <contains CONTAINER:(,,"Image Library")=SEPARATE>
    <contains CONTAINER:(,,"Image Library Group")=SEPARATE>
      <has acq context CODE:(,,"Modality")=(CT,DCM,"Computed Tomography")>
      <has acq context DATE:(,,"Study Date")="20170910">
      <has acq context UIDREF:(,,"Frame of Reference UID")="1.3.6.1.4.1.40744.29.28518703451127075549995420991770873582">

...

  <contains CONTAINER:(,,"Imaging Measurements")=SEPARATE>
    <contains CONTAINER:(,,"Measurement Group")=SEPARATE>
      <has obs context TEXT:(,,"Tracking Identifier")="Gross Target Volume">
      <has obs context UIDREF:(,,"Tracking Unique Identifier")="1.3.6.1.4.1.43046.3.1.4.0.60427.1539113881.935516"
>
      <contains CODE:(,,"Finding")=(C112913,NCIt,"Gross Target Volume")>
      <contains IMAGE:(,,"Referenced Segment")=(SG image,,1)>
      <contains UIDREF:(,,"Source series for segmentation")="1.3.6.1.4.1.40744.29.18397950185694012790332812250603
612437">
      <has concept mod CODE:(,,"Finding Site")=(T-28000,SRT,"Lung")>
      <contains NUM:(,,"shape_MeshVolume")="7.255467E+04" (1,UCUM,"no units")>
      <contains NUM:(,,"Maximum 3D diameter")="7.491328E+01" (1,UCUM,"no units")>
      <contains NUM:(,,"shape_Maximum2DDiameterSlice")="6.767570E+01" (1,UCUM,"no units")>
      <contains NUM:(,,"Elongation")="7.993260E-01" (1,UCUM,"no units")>
      <contains NUM:(,,"shape_MinorAxisLength")="4.699969E+01" (1,UCUM,"no units")>
      <contains NUM:(,,"Flatness")="6.517569E-01" (1,UCUM,"no units")>
      <contains NUM:(,,"shape_Maximum2DDiameterColumn")="6.746851E+01" (1,UCUM,"no units")>
      <contains NUM:(,,"Surface to volume ratio")="1.572168E-01" (1,UCUM,"no units")>
      <contains NUM:(,,"shape_Maximum2DDiameterRow")="6.072891E+01" (1,UCUM,"no units")>
      <contains NUM:(,,"shape_VoxelVolume")="7.285600E+04" (1,UCUM,"no units")>
      <contains NUM:(,,"Sphericity")="7.375024E-01" (1,UCUM,"no units")>
      <contains NUM:(,,"Surface area")="1.140681E+04" (1,UCUM,"no units")>
      <contains NUM:(,,"shape_MajorAxisLength")="5.879915E+01" (1,UCUM,"no units")>
      <contains NUM:(,,"shape_LeastAxisLength")="3.832275E+01" (1,UCUM,"no units")>
      <contains NUM:(,,"Small zone emphasis")="7.384502E-01" (1,UCUM,"no units")>
      <contains NUM:(,,"glszm_SmallAreaLowGrayLevelEmphasis")="3.381883E-03" (1,UCUM,"no units")>
      <contains NUM:(,,"Normalised grey level non-uniformity")="3.136554E-02" (1,UCUM,"no units")>
      <contains NUM:(,,"glszm_SmallAreaHighGrayLevelEmphasis")="5.478214E+02" (1,UCUM,"no units")>
      <contains NUM:(,,"Large zone emphasis")="3.873234E+03" (1,UCUM,"no units")>

...
```

# Questions?

* Andrey Fedorov andrey.fedorov@gmail.com

# References

* Herz C, Fillion-Robin J-C, Onken M, Riesmeier J, Lasso A, Pinter C, Fichtinger G, Pieper S, Clunie D, Kikinis R, Fedorov A. dcmqi: An Open Source Library for Standardized Communication of Quantitative Image Analysis Results Using DICOM. Cancer Research. 2017;77(21):e87–e90 http://cancerres.aacrjournals.org/content/77/21/e87
* Fedorov A, Clunie D, Ulrich E, Bauer C, Wahle A, Brown B, Onken M, Riesmeier J, Pieper S, Kikinis R, Buatti J, Beichel RR. (2016) DICOM for quantitative imaging biomarker development: a standards based approach to sharing clinical data and structured PET/CT analysis results in head and neck cancer research. PeerJ 4:e2057 https://doi.org/10.7717/peerj.2057
