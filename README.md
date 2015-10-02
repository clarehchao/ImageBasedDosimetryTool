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


### Software System Setup

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


### Toolkit Functions

Before running the following commands, be sure to do the following:

1. Set up the patient data directory
  - PT0001/Vol
  - PT0001/GeoIMdata
2. Create a .json for the patient case


###### Segment the patient CT images into user-defined organ and tumor contoured manually or automatically
```
./CT2G4files.py inputfile/________.json
```

2. Convert the segmented CT image into Geant4 input files for [Monte Carlo dosimetry evaluation](https://github.com/clarehchao/VoxelizedHumanDoseMultiSDv1) 
3. Process the Monte-Carlo simualtion output files to compute the dose factors for a given patient (S-values and etc.)
4. For each source organ, quantify the time activy curve from the PET images and estimate the residence time




