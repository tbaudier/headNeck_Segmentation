#!/usr/bin/env python

import os
import itk
import click
import gatetools as gt
import re
import subprocess
import pydicom
import json
import numpy as np
import jellyfish
import sys
sys.path.append('/home/tbaudier/Software/syd_algo/')
from anonymize import *
sys.path.append('/home/tbaudier/clbwiki/code/patientIdEncryption/')
from encryptId import *

mainDict = {
    "CTV": ["CTVNBR", "CTVNHR", "CTVTHR", "CTVTBR", "CTVPROSTATE", "CTVVS", "CTV-VS", "CTV", "CTVPROSTATEVS", "CTVT", "CTVPROFIT", "CTVN", "CTVNPELVIS", "CTVNOBT", "CTVND", "CTVNG", "CTVNLOMBO"],
    "PTV": ["PTVBRART", "PTVHRART", "PTVIRART", "PTVIR", "PTVHR", "PTVPROSTATE", "PTV78", "PTVBIS", "PTVGY", "PTVVS", "PTV", "PTVVS", "PTVPROSTATEVS", "PTVBR", "PTVPROFIT", "PTVN", "PTV66", "PTVNODULE", "PTV54"],
    "Bladder_ext": ["VESSIE", "BLADDER", "VESSIEPROFIT", "BLADDERPROFIT", "VESSIEEXT", "BLADDEREXT", "VESSIEENTIER", "BLADDERENTIER", "VESSIEEXTPROFIT", "BLADDEREXTPROFIT"],
    "Bladder_int": ["VESSIEINT", "BLADDERINT", "VESSIEINTPROFIT", "BLADDERINTPROFIT", "VESSINT"],
    "Bladder_wall": ["PAROIVESICALE", "BLADDERWALL", "PV", "PVPROFIT", "PAROIVESSIE"],
    "Rectum_ext": ["RECTUM", "RECTUMPROFIT", "RECTUMEXT", "RECTUMENTIER", "RECTUMEXTPROFIT"],
    "Rectum_int": ["RECTUMINT", "RECTUMINTPROFIT", "RECTINT"],
    "Rectum_wall": ["PAROIRECTALE", "PR", "RECTUMWALL", "PRPROFIT", "PAROIRECT"],
    "Patient": ["PATIENT", "EXTERNAL", "BODY"],
    "Prostate": ["PROSTATE"],
    "GTV": ["GTV", "GTVT", "GTVN", "GTVN-Il", "GTVTIRM1"],
    "Femoral_heads": ["TETEFEM", "FEMORALHEAD", "FEMUR", "TETEFEMORALE", "TETEFEMORALED", "TETEFEMORALEG"],
    "Kidneys": ["REIN", "KIDNEY", "REINDT", "REING"],
    "Spinal_cord": ["MOELLEEPINIERE", "MOELLE", "SPINAL"],
    "Liver": ["LIVER", "FOIE"],
    "Gastrointestinal_tract": ["INTESTINS", "BOWEL", "PERITONEAL", "DIGESTIF", "SIGMOIDE", "GRELE", "INTESTIN"],
    "OAR": ["OAR", "PRVTC3", "PRVMOELLE5", "PRV", "PRVPLEXUS3"],
    "Testicles": ["TESTICULES", "SCROTUM"],
    "Penis": ["PENIS", "VERGE"],
    "Perineum": ["BULBE", "PELVIEN", "PENILBULB"],
    "Noeud_pelvis": ["NPELVIEN", "PELVIS", "PUBIS"],
    "Seminal_vesicles": ["VESICULESSEMINALES", "VS", "VS_IRM", "VSPROX"],
    "Lymph_node": ["GANGLIONLINFATIQUE", "GG", "GANGLION", "GGPTTPELVIS", "GGPELVIENS", "GANGLIONSPELV", "GGOBTUR", "AIRESGG"],
    "Canal_anal": ["CANALANAL"],
    "Iliac": ["ILIAQUE", "IL", "ILPRIMITIF", "ILIAQINTERNE", "ILIAQEXTERNE"],
    "Undefined": ["ARTIFACT", "OINTG", "ATM", "TRUC", "BOLUS", "BILLE", "CAI", "R", "TEMP", "TRANSMITTER", "CLOU", "NAIL", "OGE", "RPTRANSMITTER", "RPCENTERTRANSMITTER", "BOITE", "PROT+OS", "LP.VB", "ADP", "BOOSTHYPOGAU", "BOOST", "INTERSECTION", "RAYPILOT", "N0", "PROTHESE", "PR.INTERNE", "SACDIG", "GRAIN", "GRAINMOY", "COTYLEG", "QUEUECHEVAL", "TF"],
    "Buccal_cavity": ["CAVITEBUCCALE"],
    "Chiasma": ["CHIASMA"],
    "Crystalline_lens": ["CRISTALLIND", "CRISTALLING"],
    "Larynx": ["LARYNX"],
    "Mandible": ["MANDIBULE", "MANDIBLE"],
    "Optic_nerve": ["NERFOPTIQUED", "NERFOPTIQUEG", "NOD", "NOG", "NOPTD", "NOPTG"],
    "Eye": ["OEILD", "OEILG", "EYEGLOBERIGHT", "EYEGLOBELEFT"],
    "Parotid_gland": ["PAROTIDEG", "PAROTIDED", "PAROTIDGLANDLEFT", "PAROTIDGLANDRIGHT"],
    "Trachea": ["TRACHEE"],
    "Brainstem": ["TRONC", "TRONCCEREBRAL"],
    "Thyroid": ["THYROIDE"],
    "Inner_ear": ["COCHLEED", "COCHLEEG", "OREILLEINTERNEG", "OREILLEINTERNED"],
    "Brachial_plexus": ["PLEXUSBRACHIALE"],
    "Esophagus": ["OESOPHAGE"],
    "Hypophysis": ["HYPOPHYSE"],
    "Lips": ["LEVRES", "LIPS"],
    "Medulllary_cavity": ["CANALMEDULLAIRE"],
    "Brain": ["ENCEPHALE", "BRAIN", "CERVEAU"]
}

def invertTransfoMat(transfoMat):
    """Invert an affine transformation matrix.
    """
    transfoMatInv = np.zeros_like(transfoMat)
    A = transfoMat[:3, :3] # the rotation part in the affine matrix
    b = transfoMat[:3,3]  # the translation part in the affine matrix
    transfoMatInv[:3,:3] = A.T # A is a rotation matrix so its transpose is its inverse
    transfoMatInv[:3,3] = -A.T.dot(b)
    transfoMatInv[3,3] = 1.
    return transfoMatInv


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split('(\d+)', text)]

def natural_keys_second(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split('(\d+)', text[1])]

def write(data, filename):
    with open(filename, 'w') as output:
        json.dump(data, output)

def load(filename):
    with open(filename, 'r') as input:
        return(json.load(input))

def distanceIdentifyStruct(roi):
    distances = {}
    minKey = list(mainDict.keys())[0]
    for key in mainDict:
        distances[key] = 1e12
        for roiModel in mainDict[key]:
            distances[key] = min(distances[key], jellyfish.levenshtein_distance(roiModel, roi))
        if distances[minKey] > distances[key]:
            minKey = key
    return(minKey)

def identifyStruct(roi):
    #return roi
    if "LOGE" in roi:
        print("Warning LOGE roi:" + roi)
    for key in mainDict:
        if roi in mainDict[key]:
            return key

    print("Warning roi:" + roi)
    minKey = distanceIdentifyStruct(roi)
    print("    Identify as :" + minKey)
    return minKey

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--inputfolder', default='.', help='Input folder where patient is present')
@click.option('--ipp', default=None, help='Id IPP research')
@click.option('--auto', is_flag=True, help='To do all registered patient automatically')
@click.option('--writejson', is_flag=True, help='Redo completly json file without the output folders')
def convertPatient(inputfolder, ipp, auto, writejson):
    """
    \b
    For inputfolder patient_xxxxxxx and ipp recherche yyyyyyyy
    If ipp research is None, create IPP with Knuth algorithm
    Auto look for all patient folder
    Create the output folder patient.yyyyyyyy in the same directory than patient_xxxxxxx containing:
        - input
          - CT.nii
          - all binary cropped structures .nii
          - anonymised rtstruct.dcm
          - anonymised plan.dcm
          - cbct.nii reconstruction for all available cbct (cbct.0.nii is the first cbct for that treatment) (if 2 cbct are at the same day, I take the second and remove the first because the first is mis-aligned)
        - cbctCreation
          - CT.nii (CT.nii aligned in the coordinate system of cbct.0.nii)
          - geom.xml
          - Gate ou rtkforwardprojection output
          - reconstruction
          - output in CT.nii coordinate system (to be able to use structures)
    """
    if not auto:
        convertPatientAPI(inputfolder, ipp, writejson)
    else:
        patientInputDirectory = "/export2/home/tbaudier/studies/cbctHeadNeck/images/train/synergy/"
        patients = []
        ipps = []
        if writejson:
          os.remove("delpel.json")
        for dir in os.listdir(patientInputDirectory):
            if dir.startswith("patient_"):
                patients += [os.path.join(patientInputDirectory, dir)]
                id = int(dir[8:])
                ipps += [str(encryptId(id))]

        for patient, ippNumber in zip(patients, ipps):
            print("Processing: " + patient + " " + ippNumber)
            convertPatientAPI(patient, ippNumber, writejson)

def convertPatientAPI(inputfolder, ipp, writejson=False):
    #Create output directory based on ipp research
    patientParentDirectory = os.path.dirname(inputfolder)
    id = inputfolder[-7:]
    outputDirectory = os.path.join(patientParentDirectory, "output", "patient." + str(ipp))
    newPatientOutputFolder = True
    if os.path.isdir(outputDirectory):
      newPatientOutputFolder = False
    if not writejson and not newPatientOutputFolder:
      print(outputDirectory + " already exists")
      return(0)
    if newPatientOutputFolder:
      os.makedirs(outputDirectory)
      os.makedirs(os.path.join(outputDirectory, "input"))

    #Find all studies irradiation in the patient folder
    studies = {}
    inputRTDose = []
    inputRTStruct = []
    for root, dirs, files in os.walk(inputfolder):
        for dir in dirs:
            if root.endswith("CT_SET"):
                studies[dir] = {}
                studies[dir]["inputCT"] = []
                studies[dir]["inputCBCT"] = []
                studies[dir]["inputRTStruct"] = ""
                studies[dir]["inputPlan"] = ""
                studies[dir]["inputINI"] = ""
                if newPatientOutputFolder:
                  os.makedirs(os.path.join(outputDirectory, "input", dir))

    #Find all input images
    for root, dirs, files in os.walk(inputfolder):
        for file in files:
            if file.startswith("CT_IMAGE_") and file.endswith(".DCM"):
                dir = os.path.basename(root)
                studies[dir]["inputCT"] += [os.path.join(root, file)]
            if file.startswith("DCMTPS_Calculated") and file.endswith(".dcm"):
                inputRTDose += [os.path.join(root, file)]
            if root.endswith("CT_SET") and file.endswith(".DCM"):
                inputRTStruct += [os.path.join(root, file)]
            if root.endswith("DICOM_PLAN") and file.endswith(".DCM"):
                ds = pydicom.read_file(os.path.join(root, file))
                dir = ds[(0x0008, 0x0018)].value
                if dir in studies.keys():
                    studies[dir]["inputPlan"] = os.path.join(root, file)
                    struct = ds[(0x300c,0x0060)][0][(0x0008,0x1155)].value
                    studies[dir]["inputRTStruct"] = struct
            if root.endswith("Reconstruction") and file.endswith(".SCAN"):
                year = file.split(".")[-2][:4]
                month = file.split(".")[-3]
                day = file.split(".")[-4]
                hour = file.split(".")[-2][4:6]
                minute = file.split(".")[-2][6:8]
                second = file.split(".")[-2][8:10]
                if os.path.isfile(os.path.join(root, file[:-4] + "INI")):
                    bashCommand = "cat " + os.path.join(root, file[:-4] + "INI")
                    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                    output, error = process.communicate()
                    output = output.split(b'\r\n')
                    for element in output:
                      if element.startswith(b"ReferenceUID="):
                        dir = element.split(b'=')[1].decode("utf-8")
                    studies[dir]["inputCBCT"] += [(os.path.join(root, file), str(year)+str(month)+str(day)+str(hour)+str(minute)+str(second))]
            if "CT_SET" in root and file.startswith("1.") and file.endswith(".INI"):
                dir = os.path.basename(root)
                studies[dir]["inputINI"] = os.path.join(root, file)

    for struct in inputRTStruct:
        ds = pydicom.read_file(struct)
        structId = ds[(0x0008, 0x0018)].value
        for study in studies.keys():
            if structId == studies[study]["inputRTStruct"]:
                studies[study]["inputRTStruct"] = struct

    for dir in studies.keys():
        studies[dir]["inputCT"].sort(key=natural_keys)
        studies[dir]["inputCBCT"].sort(key=natural_keys_second)

        #Remove wrong cbct (if 2 cbct are done the same day, keep the second -> it means that the first has a wrong alignment)
        removeIndexCBCT = []
        for indexCBCT in range(1, len(studies[dir]["inputCBCT"])):
          if int(studies[dir]["inputCBCT"][indexCBCT][1]) - int(studies[dir]["inputCBCT"][indexCBCT-1][1]) < 10000:
            removeIndexCBCT += [indexCBCT-1]
        studies[dir]["inputCBCT"] = [i for j, i in enumerate(studies[dir]["inputCBCT"]) if j not in removeIndexCBCT]

        #Convert CT
        inputDicomCT = gt.read_dicom(studies[dir]["inputCT"])
        if newPatientOutputFolder:
          #outputCT = gt.image_convert(inputDicomCT)
          itk.imwrite(inputDicomCT, os.path.join(outputDirectory, "input", dir, "CT.nii"))
        transfoCT = np.array([[1.,0,0,0], [0,1.,0,0], [0,0,1.,0], [0,0,0,1.]])

        #Convert CBCT threw clitk
        cbctIndex = 0
        for cbct in studies[dir]["inputCBCT"]:
          if os.path.isfile(cbct[0][:-4] + "INI.XVI"):
            #Read CBCT transformation matrix
            with open(cbct[0][:-4] + "INI.XVI") as f:
              for line in f.readlines():
                  if "OnlineToRefTransformCorrection=" in line:
                      outputCBCT = os.path.join(outputDirectory, "input", dir, "cbct." + str(cbctIndex) + ".nii")
                      tempLine = line.split('=')[1]
                      if tempLine[0] == " ":
                          tempLine = tempLine[1:]
                      transfoMatrix = np.array([float(t) for t in tempLine.split()]).reshape(4,4).T
                      transfoMatrix2 = np.copy(transfoMatrix)
                      transfoMatrix2[:-1, -1] *= 10 #convert from cm to mm
                      transfoMatCT2CBCT = transfoCT.dot(transfoMatrix2)
                      transfoMatCBCT2CT = invertTransfoMat(transfoMatCT2CBCT)
                      np.savetxt(outputCBCT[:-3] + "mat", transfoMatCBCT2CT)
                      if newPatientOutputFolder:
                        bashCommand = "clitkImageConvert -i " + cbct[0] + " -o " + outputCBCT
                        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                        output, error = process.communicate()
                        bashCommand = "clitkAffineTransform -i " + outputCBCT + " -o " + outputCBCT + " -m " + outputCBCT[:-3] + "mat" + " --pad=-1024 --transform_grid"
                        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                        output, error = process.communicate()
                      os.remove(outputCBCT[:-3] + "mat")
                      np.savetxt(outputCBCT[:-3] + "mat", transfoMatrix)
                      break
            cbctIndex += 1



        #Copy Frame.xml
        cbctIndex = 0
        for cbct in studies[dir]["inputCBCT"]:
            files = []
            for file in  os.listdir(os.path.join(os.path.dirname(cbct[0]), "..")):
                if os.path.isfile(os.path.join(os.path.dirname(cbct[0]), "..", file)) and "Frames" in file:
                    files += [os.path.join(os.path.dirname(cbct[0]), "..", file)]
            if os.path.isdir(os.path.join(os.path.dirname(cbct[0]), "..", "..", "..", "ARCHIVE")):
                for file in  os.listdir(os.path.join(os.path.dirname(cbct[0]), "..", "..", "..", "ARCHIVE")):
                    if os.path.isfile(os.path.join(os.path.dirname(cbct[0]), "..", "..", "..", "ARCHIVE", file)) and file.startswith(os.path.basename(os.path.dirname(os.path.dirname(cbct[0])))[4:]):
                        files += [os.path.join(os.path.dirname(cbct[0]), "..", "..", "..", "ARCHIVE", file)]

            if len(files) == 0:
                print("NO FRAMES FOR: " + os.path.basename(os.path.dirname(os.path.dirname(cbct[0])))[4:])
            
            #copy the frames in the output folder and modify the patient information
            fileIndex = 0
            for file in files:
                fold = open(file, "r")
                fnew = open(file + ".tmp", "w")
                isPatient = False
                for line in fold:
                    if "<Frames />" in line:
                        print("Wrong patient frame")
                    if "<Patient>" in line:
                        isPatient = True
                        fnew.write(line)
                        fnew.write("    <FirstName> </FirstName>\n")
                        fnew.write("    <MiddleName> </MiddleName>\n")
                        fnew.write("    <LastName> </LastName>\n")
                        fnew.write("    <ID>" + str(ipp) + "</ID>\n")
                    elif "</Patient>" in line:
                        isPatient = False
                        fnew.write(line)
                    elif isPatient:
                        continue
                    else:
                        fnew.write(line)
                fold.close()
                fnew.close()
                shutil.move(file + ".tmp", os.path.join(outputDirectory, "input", dir, "cbct." + str(cbctIndex) + "_" + str(fileIndex) + ".xml"))
                fileIndex += 1
            cbctIndex += 1

        #Convert struct threw clitk and after to .nii
        inputTmpRTStruct = []
        if newPatientOutputFolder and not studies[dir]["inputRTStruct"] == "":
          bashCommand = "clitkDicomRTStruct2Image -c --mha -i " + studies[dir]["inputRTStruct"] + " -o " + os.path.join(outputDirectory, "input", dir, "tmp.rtstruct.") + " -j " + os.path.join(outputDirectory, "input", dir, "CT.nii")
          process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
          output, error = process.communicate()
          for root, dirs, files in os.walk(os.path.join(outputDirectory, "input", dir)):
            for file in files:
              if file.startswith("tmp.rtstruct.") and file.endswith(".mha"):
                if not file == file.encode('utf-8', 'replace').decode():
                    continue
                inputTmpRTStruct += [file]

        structName = []
        if newPatientOutputFolder:
          for struct in inputTmpRTStruct:
            inputStruct = itk.imread(os.path.join(outputDirectory, "input", dir, struct))
            itk.imwrite(inputStruct, os.path.join(outputDirectory, "input", dir, struct[13:-3] + "nii"))
            structName += [struct[13:-4]]
            os.remove(os.path.join(outputDirectory, "input", dir, struct))
        else:
          for root, dirs, files in os.walk(os.path.join(outputDirectory, "input", dir)):
            for file in files:
              if not file.startswith("cbct.") and not file.startswith("CT.") and file.endswith(".nii"):
                structName += [file[:-4]]
          

        #Anonymise rtstruct
        if newPatientOutputFolder and not studies[dir]["inputRTStruct"] == "":
            anonymizeDicomFile(studies[dir]["inputRTStruct"], os.path.join(outputDirectory, "input", dir, "rtstruct.dcm"), "anonymous", ipp, False, [])

        #Anonymise dicom plan
        if newPatientOutputFolder and not studies[dir]["inputPlan"] == "":
          anonymizeDicomFile(studies[dir]["inputPlan"], os.path.join(outputDirectory, "input", dir, "rtplan.dcm"), "anonymous", ipp, False, [])

        #Anonymise rt dose
        if newPatientOutputFolder:
          indexDose = 0
          for dose in inputRTDose:
            anonymizeDicomFile(dose, os.path.join(outputDirectory, "input", dir, "rtdose" + str(indexDose) + ".dcm"), "anonymous", ipp, False, [])
            indexDose += 1

        #Keep treatment from .INI file
        bashCommand = "cat " + studies[dir]["inputINI"]
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        output = output.split(b'\r\n')
        for element in output:
          if element.startswith(b"TreatmentID="):
            treatmentName = element.split(b'=')[1].decode("utf-8")

        #Write id and ipp in json file
        try:
            jsonData = load('delpel.json')
        except:
            jsonData = {}
            jsonData['patients'] = {}
            jsonData['header'] = {
                "key": "ipp research",
                "treatment id": "TreatmentID value found in .INI file",
                "clarity": "Does the patient treated with Clarity (artifacts)",
                "machine": "Name of the treatment machine",
                "city": "Treatment city",
                "ROI name in rtStruct": "found corresponding name"
            }
        if not ipp in jsonData['patients']:
            jsonData['patients'][ipp] = {}
        if not dir in jsonData['patients'][ipp]:
          jsonData['patients'][ipp][dir] = {
                              'treatment id': treatmentName,
                              'ReferenceUID': dir,
                              'Structures': {},
                              'clarity': False,
                              'city': "Lyon",
                              'machine': os.path.basename(os.path.dirname(inputfolder))
                              }
          for struct in structName:
              jsonData['patients'][ipp][dir]['Structures'][struct] = identifyStruct("_".join(struct.split("_")[1:]).upper())

          with open("delpel.json", "w") as jsonFile:
              json.dump(jsonData, jsonFile, indent=4, sort_keys=True)

if __name__ == '__main__':
    convertPatient()
