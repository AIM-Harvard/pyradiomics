from __future__ import annotations

import csv
import os

SqDic = {}
SqDic["T2-TSE-TRA"] = "t2tra"
SqDic["T2-TRA"] = "t2tra"
SqDic["T2-TSE-SAG"] = "t2sag"
SqDic["T2-SAG"] = "t2sag"
SqDic["T2-TSE-COR"] = "t2cor"
SqDic["T2-COR"] = "t2cor"
SqDic["B1000"] = "dwi"
SqDic["B1100"] = "dwi"
SqDic["ADC"] = "adc"

LabelDic = {}
# LabelDic['Reader-1'] = '2 Semi-Auto_Inexp-3'
# LabelDic['Reader-2'] = '2 Semi-Auto_Inexp-4'
LabelDic["Reader-3"] = "3 Semi-Auto_Exp-1"
LabelDic["Reader-4"] = "3 Semi-Auto_Exp-2"
# LabelDic['Reader-5'] = '4 Manual-1'
# LabelDic['Reader-6'] = '4 Manual-2'
LabelDic["Reader-7"] = "1 Auto-1"
LabelDic["Reader-8"] = "3 Semi-Auto_Exp-1"
LabelDic["Reader-9"] = "2 Semi-Auto_Inexp-1"
LabelDic["Reader-10"] = "2 Semi-Auto_Inexp-2"
LabelDic["Reader-11"] = "3 Semi-Auto_Exp-2"


def main():
    DATA_ROOT_PATH = r"T:/Research/07. Current Projects/2. Robust Radiomics/1. Slicer Dataset COMPLETE/"
    inputDirectory = DATA_ROOT_PATH + r"/Included"
    outputFile = DATA_ROOT_PATH + r"/Included/FileList.csv"
    filetype = ".nrrd"

    print("Scanning files...")

    datasetHierarchyDict = scanpatients(inputDirectory, filetype)

    print(f"Found {len(datasetHierarchyDict.keys())} patients, writing csv")

    try:
        with open(outputFile, "w") as outFile:
            cw = csv.writer(outFile, lineterminator="\n")
            cw.writerow(["Patient", "Sequence", "Reader", "Image", "Mask"])

            for patient, Studies in sorted(
                datasetHierarchyDict.items(), key=lambda t: t[0]
            ):
                for Study, im_fileList in sorted(
                    Studies["reconstructions"].items(), key=lambda t: t[0]
                ):
                    for i_idx, im_file in enumerate(im_fileList):

                        if Study in Studies["segmentations"]:
                            for Reader, seg_fileList in sorted(
                                Studies["segmentations"][Study].items(),
                                key=lambda t: t[0],
                            ):
                                for s_idx, seg_file in enumerate(sorted(seg_fileList)):

                                    i_name = Study
                                    if i_idx > 0:
                                        i_name += f" ({i_idx + 1!s})"

                                    s_name = Reader
                                    if s_idx > 0:
                                        s_name += f" ({s_idx + 1!s})"

                                    cw.writerow(
                                        [patient, i_name, s_name, im_file, seg_file]
                                    )
    except Exception as exc:
        print(exc)


def scanpatients(f, filetype):
    outputDict = {}

    for dirpath, _dirnames, filenames in os.walk(f):
        # list of all filenames. check pt number, check if roi, check sq
        for fname in filenames:
            if (fname[0:3] == "Pt-") & (fname.endswith(filetype)):
                PtNo = fname[3:7]

                if PtNo not in outputDict:
                    outputDict[PtNo] = {"reconstructions": {}}
                    outputDict[PtNo]["segmentations"] = {}

                for SqKey, SqVal in SqDic.items():
                    if ("ROI_" + SqVal) in fname:
                        for ReaderKey, ReaderVal in LabelDic.items():
                            if (ReaderKey + "_") in fname:
                                if SqVal not in outputDict[PtNo]["segmentations"]:
                                    outputDict[PtNo]["segmentations"][SqVal] = {}
                                if (
                                    ReaderVal
                                    not in outputDict[PtNo]["segmentations"][SqVal]
                                ):
                                    outputDict[PtNo]["segmentations"][SqVal][
                                        ReaderVal
                                    ] = set()
                                outputDict[PtNo]["segmentations"][SqVal][ReaderVal].add(
                                    os.path.join(dirpath, fname)
                                )
                                break
                    elif SqKey in fname:
                        if SqVal not in outputDict[PtNo]["reconstructions"]:
                            outputDict[PtNo]["reconstructions"][SqVal] = set()
                        outputDict[PtNo]["reconstructions"][SqVal].add(
                            os.path.join(dirpath, fname)
                        )
                        break
    return outputDict


if __name__ == "__main__":
    main()
