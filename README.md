## ImageBasedDosimetryTool
A Python Software toolkit to evaluate image-based dosimetry for targeted radionuclide therapy

### Some key words to define:

- Source organ: the organ where radioactivity originates from (SPECT/PET nuclear imaging: the image regions that uptake the injected radiopharmaceutical)
- Target organ: the organ that the radiation/radioactivity from the source organ deposits further radiation dose to
- S-value: the dose factor estimated for a given patient geometry from a source organ to a target organ
- Time activity curve: a curve that relates the time since the initial tracer injection with % of injected dose
- Residence time: the duration of radiation in a given organ, i.e. total aread under the time activity curve

### Reference for nuclear medicine dosimetry:

[Fundamentals of Nuclear Medicine Dosimetry by Michael G. Stabin] (http://www.springer.com/us/book/9780387745787)


Software System Setup
--------------------

- Python 2.7.10, Anaconda 2.3.0 (64-bit) for Mac OS or Ubuntu Linux
- on UCSF PRL Higgs server, in your personal account, make sure to add the following line in the ~/.bashrc file
```
# added by Anaconda 2.3.0 installer
export PATH="/data1/anaconda/bin:$PATH"
```
- the toolkit was built in the Python 2.7 environment, so be sure to invoke the py27 environment before running the code:
```
# Start the python 2.7 environment
source activate py27

# Exit the python 2.7 environment
source deactivate
``` 

If any case, you are not using PRL Higgs server to run this toolkit, one can still run this toolkit by [installing Anaconda python] (https://www.continuum.io/downloads) and be sure to install the following packages via binstar or conda install:
- pydicom
- mysql-python
- lmfit


Toolkit Functions
--------------------

Before running the following commands, be sure to do the following:

1. Set up the patient data directory [PTdir], e.g. pt_id = 1
  - PT0001/GeoIMdata/[geo_id]: the raw pixel files of all the contoured organs and organs (.tsv) or .bin (if a binary files for a specific contour is already created)
  - PT0001/VOIs_Amide: Amide contour files (.xif) for all imaging acquisitions and the ROI statistics of all the contoured VOIs (.tsv)
  - PT0001/PMODData (optional): If one were to use PMOD to contour, this directory includes all the PMOD VOI statistics files (.voistat) for eachh organ at all imaging time points (e.g. PT0001/PMODData/Tumor1/_____.voistat)
  - PT0001/IM: dicom images of all the PET/CT images acquired at several time points
  - PT0001/Summary: the output dose plot, residence time and mass vs. target organ plot, and PET/CT images of the patient with identified tumors

2. Create a .json for the patient case
  - "G4dir": Geant4 data directory, e.g. "/data2/G4data_Clare",
  - "VHDMSDdir": toolkit code directory, e.g. "/data1/Proj_I131MIBGTherapy/VHDMSDlite"
  - "NBptdir": patient directory, e.g. "/data1/Proj_I131MIBGTherapy/ptData",
  - "fwdir": Geant4 input file WRITE directory, e.g. "/data2/G4data_Clare/G4INPUT",
  - "ctdir": the dicom image directory of the CT images that the organ contour was based on for CT segmentation, e.g. "/data1/Proj_I131MIBGTherapy/ptData/PT0002/IM/E13450/4",
  - "geotag": Geant4 geoemtry name, e.g. "segCT_MIBGPT2",
  - "PTid": the patient ID
  - "pmodftag": the format of PMOD .voistat files
  - "img_isotope_halflife_day": the half life of the imaging isotope (unit: days)
  - "therapy_isotope_halflife_day":the half-life of the therapy isotope (unit: days)
  - "I131MIBGinjDosemCi": the injected dose of the I-131 MIBG therapy (unit: mCi) 
  - "simpkg": Geant4 version, e.g. "G4.9.6.p02"
  - "isDosePlot": plot and save the dose-related figures if true, otherwise
  - "isFitPlot": plot the fitted time activity curves if true, otherwise
  - "G4AppName": Geant4 application name
  - "G4SimDataDir": Geant4 simulation data directory and computer name
  - "frtype": binary image data type for file read, e.g."uint8"
  - "frindx":[2,3,4]
  - "ecomptag":"adf"
  - "phantomtag":""
  - "masktag":[]
  - "binfwtype":"uint8"
  - "nxyz":[512,512,364]
  - "dxyz":[1.36719,1.36719,5.0]
  - "xyz0":[-357.10,-315.70,-30.50]
  - "HUthresh":[-5000,-400,200,1440,5000]
  - "HUthresh_name":["Air(inbody)","ResidualSoftTissue","Cranium","Teeth"]
  - "organvoiname":["Lung","Brain","Heart","Liver","SalivaryGlands","Spleen","Urinarybladder","Stomach","Thyroid","Kidney"]
  - "tumorvoiname":["Tumor1","Tumor2","Tumor3","Tumor4","Tumor5","Tumor6","Tumor7","Tumor8","Tumor9"]
  - "srcparticle":"I131"
  - "excludesrcname":[]
  - "therun":[1,10,1]
  - "srcname":{"TotalBody":["ResidualSoftTissue","Cranium","Teeth","Lung","Brain","Heart","Liver","SalivaryGlands","Spleen","Urinarybladder","Stomach","Thyroid","Kidney","Tumor1","Tumor2","Tumor3","Tumor4","Tumor5","Tumor6","Tumor7","Tumor8","Tumor9"],"Heart":["Heart"],"Liver":["Liver"],"SalivaryGlands":["SalivaryGlands"],"Spleen":["Spleen"],"Urinarybladder":["Urinarybladder"],"Stomach":["Stomach"],"Thyroid":["Thyroid"],"Kidney":["Kidney"],"Tumor1":["Tumor1"],"Tumor2":["Tumor2"],"Tumor3":["Tumor3"],"Tumor4":["Tumor4"],"Tumor5":["Tumor5"],"Tumor6":["Tumor6"],"Tumor7":["Tumor7"],"Tumor8":["Tumor8"],"Tumor9":["Tumor9"],"Brain":["Brain"],"Lung":["Lung"]}}


##### Segment & Convert CT images to Geant4 input files

- Segment the patient CT images into user-defined organ and tumor contoured manually or automatically
- Convert the segmented CT image into Geant4 input files for [Monte Carlo dosimetry evaluation](https://github.com/clarehchao/VoxelizedHumanDoseMultiSDv1) 
```
./CT2G4files.py inputfile/________.json
```
Note: if [G4iniputdir]/GeometryIM/binIM/[geo_id]/GeoVol.bin does exist, the code will only create the source map defined in 'srcname' in the .json file.

##### Compute the Dose Factors from Monte Calor Simulations
- Process the Monte-Carlo simualtion output files to compute the dose factors for a given patient (S-values and etc.)
- Save the S-value from all source-to-target organ pairs into the MySQL database (UCSFDoseDB: DoseInfo table)
```
./getSvalue_mysql.py inputfile/________.json
```

##### Compute the Organ Mass and Volume of a given patient geometry
- Compute the mass and volume of all source organs defined in .json
- Save the organ mass, volume and name into the MySQL database (UCSFDoseDB: GeoInfo table)
```
./OrganMass_mysql.py inputfile/________.json
```

#### PET-image Residence Time Evaluation
- For each source organ, quantify the time activy curve from the PET images and estimate the residence time
- Save the residence time data (bi-exponential fit parameters and residence time) into the MySQL database (UCSFDoseDB: ResidenceTimeInfo table)
- Save the patient data info (pt_id, therapy dose, etc.) in the MySQL database (UCSFDoseDB: MIBGPTInfo table)
- Save the total absorbed dose to each organ/tissue in the MySQL database (UCSFDoseDB: AbsorbedDoseInfo table)
- Plot the final absorbed dose, residence time, and organ mass for each target organ and save in the directory [PTdir]/Summary/
```
./ResTime_mysql.py inputfile/________.json
```




