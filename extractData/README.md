
**Extract data from xvi**

1/ findHeandNeck.sh:  
Get the list of patient with the data and a treatment name "arynx" for larynx or pharynx

2/ createPatient.py:  
Convert the original patient to anonymized patient. Convert the dicom to .nii image, get the geom.xml and structure.
Try to guess the name of the structure

3/ correspondance.json:
Output of createPatient.py with the patient guessed stuctures
 
