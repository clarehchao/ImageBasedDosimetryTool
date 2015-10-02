## ImageBasedDosimetryTool
A Python Software toolkit to evaluate image-based dosimetry for targeted radionuclide therapy

##### The toolkit performs the following functions:

1. Segment the patient CT images into user-defined organ and tumor contoured manually or automatically
2. Convert the segmented CT image into Geant4 input files for [Monte Carlo dosimetry evaluation](https://github.com/clarehchao/VoxelizedHumanDoseMultiSDv1) 
3. Process the Monte-Carlo simualtion output files to compute the dose factors for a given patient (S-values and etc.)
4. For each source organ, quantify the time activy curve from the PET images and estimate the residence time

##### Some key words to define:

- Source organ: the organ where radioactivity originates from (SPECT/PET nuclear imaging: the image regions that uptake the injected radiopharmaceutical)
- Target organ: the organ that the radiation/radioactivity from the source organ deposits further radiation dose to
- S-value: the dose factor estimated for a given patient geometry from a source organ to a target organ
- Time activity curve: a curve that relates the time since the initial tracer injection with % of injected dose
- Residence time: the duration of radiation in a given organ, i.e. total aread under the time activity curve

##### Reference for nuclear medicine dosimetry: [Fundamentals of Nuclear Medicine Dosimetry by Michael G. Stabin] (http://www.springer.com/us/book/9780387745787)


### Software System Setup



