## ImageBasedDosimetryTool
A Python Software toolkit to evaluate image-based dosimetry for targeted radionuclide therapy

### Some key words to define:

- Source organ: the organ where radioactivity originates from (SPECT/PET nuclear imaging: the image regions that uptake the injected radiopharmaceutical)
- Target organ: the organ that the radiation/radioactivity from the source organ deposits further radiation dose to
- S-value: the dose factor estimated for a given patient geometry from a source organ to a target organ
- Time activity curve: a curve that relates the time since the initial tracer injection with % of injected dose
- Residence time: the duration of radiation in a given organ, i.e. total aread under the time activity curve

### Reference for nuclear medicine dosimetry:

[Fundamentals of Nuclear Medicine Dosimetry by Michael G. Stabin](http://www.springer.com/us/book/9780387745787)


Software System Setup
--------------------

- Python 2.7.14.final.0 Anaconda 4.4.10 (64-bit) for Mac OS or Ubuntu Linux
- On UCSF PRL Terra server, in your personal account, make sure to add the following line in the ~/.bashrc file

```
# added by Anaconda2 4.4.10 installer
export PATH="/data1/packages/anaconda2/bin:$PATH"
```
- Activate the following anaconda environment before running the code on UCSF PRL Terra server

```
source activate py27-root
```

If any case, you are not using PRL server to run this toolkit, one can still run this toolkit by [installing Anaconda python](https://www.continuum.io/downloads) and be sure to install the following packages via binstar or conda install:
- pydicom
- mysql-python
- lmfit
- pynrrd
- rootpy (installed on Linux Ubuntu 16.04 via ```conda install -c https://conda.anaconda.org/nlesc root rootpy```)


Toolkit Functions
--------------------

Before running the following commands, be sure to do the following:

1. Segment all major organs and tumors from the patient PET/CT images
  - Major organs include: brain, salivary glands, thyroid, kidney, lungs, heart, adrenal glands (if exists), spleen, liver, stomach, urinary bladder
  - Import PET and CT DICOM images into [3DSlicer](https://www.slicer.org/) and perform semi-automatic segmentation of each organ by placing seeds appropriately (see example [here](Seg_demo/3dslicer_demo_pt10.jpg))
  - Once the segmentation is complete, save and export the segmented volume into a .nrrd file
  - One would only need to segment organs from one time-point of the PET/CT imaging
  - In order import 3DSlicer-segmented volume into Amide, one need to convert the .nrrd file to DICOM images via NRRD2DCM.py
  - Once all the segmented organs are imported into Amide, compute the PET signal statistics using ROI quantification tool via Amide (Tools > Calculate ROI statistics > Execute) and save the output to .tsv (save one .tsv for each imaging time point, see an example [here](Seg_demo/PET_VOI_Quantification2.png))
  - Remember to manually place a cylinder VOI in Amide to include the entire body of the patient for PET signal quantification of the total body (see an example [here](Seg_demo/PET_VOI_Quantification1.png))
  - To get the raw pixel output of each organ from the CT image, compute the CT signal statistics using ROI quantification tool via Amide and click 'Save Raw Values' once the ROI stats is completed. See an example [here](Seg_demo/CT_VOI_RawValue.png).


2. Set up the patient data directory [PTdir], e.g. pt_id = 1
  - PT0001/GeoIMdata/[geo_id]: the raw pixel files of all the segmented organs (roi_raw_data{xxOrganName}.tsv) or .bin (if a binary files for a specific contour is already created)
  - PT0001/VOIs_Amide: Amide contour files (.xif) for all imaging acquisitions and the ROI statistics of all the contoured VOIs (.tsv)
  - PT0001/PMODData (optional): If one were to use PMOD to contour, this directory includes all the PMOD VOI statistics files (.voistat) for eachh organ at all imaging time points (e.g. PT0001/PMODData/Tumor1/_____.voistat)
  - PT0001/IM: dicom images of all the PET/CT images acquired at several time points
  - PT0001/Summary: the output dose plot, residence time and mass vs. target organ plot, and PET/CT images of the patient with identified tumors

3. Create a .json for the patient case
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
  - "xyz0": the origin of the 3D position of the CT image volume in Amide, unit: mm (can be calculated by xyz_center - dxyz*(nxyz-1)*0.5, where xyz_center is found in the 'Center' Tag when right-click on the CT image series in Amide)
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
- Convert the segmented CT image into Geant4 input files for [Monte Carlo dosimetry evaluation](https://git.radiology.ucsf.edu/PRL/VoxelizedMonteCarloDosimetry)

```
./CT2G4files.py inputfile/________.json
```
Note: if [G4iniputdir]/GeometryIM/binIM/[geo_id]/GeoVol.bin does exist, the code will only create the source map defined in 'srcname' in the .json file.

##### Compute the Dose Factors from Monte Calor Simulations
- Process and save the S-value from all source-to-target organ pairs into the MySQL database (UCSFDoseDB: DoseInfo table)
- Run the below command after finishing running all the needed S-value Monte Carlo simulation
- No need to convert the output .root file to .dat file since rootpy python package allows direct read/access to a .root file


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

##### PET-image Residence Time Investigation
- Examine the need for evaluating residence time from 4 imaging time
- Generate scatter plot of Residence Time vs Organ among all patients
- Generate scatter plot of Absorbed Dose vs Organ among all patients
- Compute the slope of two imaging time points of PET-based Time Activity Curve (TAC)

```
./ResTime_Analysis.py [data_DIR] [DB_auth_DIR]
```
[data_DIR]: the directory where all files and figures are saved
[DB_auth_DIR]: the directory where the MySQL database authentication info is


##### MySQL Database
- All the simulation and dose data are stored in a MySQL database named 'UCSFDoseDB' on UCSF PRL Terra server
- The username and password to the database can be found in a user-based file
- To interact with the database via MySQL

```
> myql -u root -p
[type in password]

mysql> use UCSFDoseDB;
mysql> show tables; 
```
- To get all the data in a table, e.g. SimInfo,

```
select * from Siminfo; 
```
- To query a table of absorbed dose, residence time, organ mass and volume for pt_id=6 and geo_id=segCT_MIBGTPT6,

```
select c.OrganName,c.ResTime_BqhrPerBq,c.AbsorbedDose_mGy,d.Volume_cm3,d.Mass_g from (select a.pt_id,a.OrganName,ResTime_BqhrPerBq,AbsorbedDose_mGy from (ResTimeInfo a INNER JOIN AbsorbedDoseInfo b ON a.pt_id=b.pt_id and a.OrganName=b.TargetOrgan) where a.pt_id=6) as c INNER JOIN (select * from GeoInfo where geo_id='segCT_MIBGPT6') as d ON c.OrganName=d.OrganName;
```

- To save the above query into a dataframe (if using Python Pandas package)

```python
import pandas as pd
import MySQLdb as mdb
from DBTool import DoseDB as ddb

DB_auth_dir = 'xxxx'
DB_usr = 'root'
DB_pw = ddb.get_DB_auth_info(DB_auth_dir,DB_usr)
DB_name = 'UCSFDoseDB'
DB_host = '127.0.0.1'
con = mdb.connect(host=DB_host, user=DB_usr,passwd=DB_pw, db=DB_name)

qr = 'c.OrganName,c.ResTime_BqhrPerBq,c.AbsorbedDose_mGy,d.Volume_cm3,d.Mass_g from (select a.pt_id,a.OrganName,ResTime_BqhrPerBq,AbsorbedDose_mGy from (ResTimeInfo a INNER JOIN AbsorbedDoseInfo b ON a.pt_id=b.pt_id and a.OrganName=b.TargetOrgan) where a.pt_id=6) as c INNER JOIN (select * from GeoInfo where geo_id=\'segCT_MIBGPT6\') as d ON c.OrganName=d.OrganName;'
df = pd.read_sql(qr,con)
```






