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

- Python 2.7.12 Anaconda 4.2.0 (64-bit) for Mac OS or Ubuntu Linux
- on UCSF PRL Higgs server, in your personal account, make sure to add the following line in the ~/.bashrc file
```
# added by Anaconda2 4.2.0 installer
export PATH="/data1/packages/anaconda2/bin:$PATH"
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
  - "petdirs": the image study directory and series number of the PET images 
  - "geotag": Geant4 geoemtry name, e.g. "segCT_MIBGPT2",
  - "PTid": the patient ID
  - "pmodftag": the format of PMOD .voistat files
  - "img\_isotope\_halflife_day": the half life of the imaging isotope (unit: days)
  - "therapy\_isotope\_halflife_day":the half-life of the therapy isotope (unit: days)
  - "I131MIBGinjDosemCi": the injected dose of the I-131 MIBG therapy (unit: mCi) 
  - "simpkg": Geant4 version, e.g. "G4.9.6.p02"
  - "isDosePlot": plot and save the dose-related figures if true, otherwise
  - "isFitPlot": plot the fitted time activity curves if true, otherwise
  - "G4AppName": Geant4 application name
  - "G4SimDataDir": Geant4 simulation data directory and computer name
  - "frtype": binary image data type for file read, e.g."uint8"
  - "frindx": the column index to read from the Amide raw pixel files for building an binary marsk, e.g. [2,3,4]
  - "ecomptag": the phantom age usd for elemental composition defined in the Monte Carlo simulation geometry, 
    - "00f" or "00m: newborn female/male
    - "01f" or '01m: 1-year-old female/male
    - "05f" or "05m": 5-year-old female/male
    - "10f" or "10m": 10-year-old female/male
    - "15f" or 15m": 15-year-old female/male
    - "adf" or "adm": adult female/male
  - "binfwtype": file format to the binary mask volume, e.g.: "uint8"
  - "nxyz": the 3D dimension of the CT image volume, e.g. [512,512,364]
  - "dxyz": the 3D voxel size of the CT image volume, unit: mm, e.g.: [1.36719,1.36719,5.0]
  - "xyz0": the initial 3D position of the CT image volume in Amide, e.g. [-357.10,-315.70,-30.50], should be the (x,y,z) coordinate of the top left corner 
  - "HUthresh": the HU thresholds for the initial segmentation, e.g.: [-5000,-400,200,1440,5000]
  - "HUthresh_name": the name of the materials segmented using "HUthresh", e.g. ["Air(inbody)","ResidualSoftTissue","Cranium","Teeth"]
  - "organvoiname": a list of organ name, e.g.: ["Lung","Brain","Heart"]
  - "tumorvoiname": a list of tumor name, e.g.: ["Tumor1","Tumor2"]
  - "srcparticle": the source particle simulated in the Monte Carlo simulations, e.g.: "I131"
  - "excludesrcname": a list of the source name to exclude, []
  - "therun": the Monte Carlo the start Run #, the end Run #, and the Run increment for evaluating S-value's, e.g. [1,10,1]
  - "srcname": a dictionary of the source organs, the dictionary key is the name of the source organ, and the value is the organ labele associated to the source organ, e.g. {"TotalBody":["ResidualSoftTissue","Cranium","Teeth","Lung","Brain","Heart","Tumor1","Tumor2"],"Heart":["Heart"],"Tumor1":["Tumor1"],"Tumor2":["Tumor2"],"Brain":["Brain"],"Lung":["Lung"]}}


##### Segment & Convert CT images to Geant4 input files

- Segment the patient CT images into user-defined organ and tumor contoured manually or automatically
- Convert the segmented CT image into Geant4 input files for [Monte Carlo dosimetry evaluation](https://github.com/clarehchao/VoxelizedMonteCarloDosimetry) 
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

##### PET-image Residence Time Evaluation
- For each source organ, quantify the time activy curve from the PET images and estimate the residence time
- Save the residence time data (bi-exponential fit parameters and residence time) into the MySQL database (UCSFDoseDB: ResidenceTimeInfo table)
- Save the patient data info (pt_id, therapy dose, etc.) in the MySQL database (UCSFDoseDB: MIBGPTInfo table)
- Save the total absorbed dose to each organ/tissue in the MySQL database (UCSFDoseDB: AbsorbedDoseInfo table)
- Plot the final absorbed dose, residence time, and organ mass for each target organ and save in the directory [PTdir]/Summary/
```
./ResTime_mysql.py inputfile/________.json
```




