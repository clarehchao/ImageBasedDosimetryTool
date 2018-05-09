import json
import pandas as pd
from pydicom.dataset import Dataset, FileDataset
import numpy as np
import datetime, time
import nrrd
import os, errno
import copy

def LoadInputParameter(fname):
    """
    load the input file and read the parameters into a dictionary via JSON
    """
    thedict = {}
    jdec = json.JSONDecoder()
    ff = open(fname)
    for line in ff:
        thedict = jdec.decode(line)
    return thedict 
    
    
def Dict2DF(thedict,thecolname,fname):
    """
       convert a dictionary to a data frame for file saving
    """
    df = pd.DataFrame(thedict.values(),index=thedict.keys(),columns=thecolname)
    df.to_csv(fname)
    print 'Save the dataframe to {}'.format(fname)
    return df


def array2Dicom(slice_idx, pixel_array, dct_info, filename):
    """
    INPUTS:
    pixel_array: 2D numpy ndarray.  If pixel_array is larger than 2D, errors.
    filename: string name for the output file.
    """

    ## This code block was taken from the output of a MATLAB secondary
    ## capture.  I do not know what the long dotted UIDs mean, but
    ## this code works.
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = 'Secondary Capture Image Storage'
    file_meta.MediaStorageSOPInstanceUID = '1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
    file_meta.ImplementationClassUID = '1.3.6.1.4.1.9590.100.1.0.100.4.0'
    ds = FileDataset(filename, {}, file_meta=file_meta, preamble="\0" * 128)
    ds.Modality = 'CT'
    ds.ContentDate = str(datetime.date.today()).replace('-', '')
    ds.ContentTime = str(time.time())  # milliseconds since the epoch
    ds.StudyInstanceUID = '1.3.6.1.4.1.9590.100.1.1.124313977412360175234271287472804872093'
    ds.SeriesInstanceUID = '1.3.6.1.4.1.9590.100.1.1.369231118011061003403421859172643143649'
    ds.SOPInstanceUID = '1.3.6.1.4.1.9590.100.1.1.111165684411017669021768385720736873780'
    ds.SOPClassUID = 'Secondary Capture Image Storage'
    ds.SecondaryCaptureDeviceManufctur = 'Python 2.7.14'

    ## These are the necessary imaging components of the FileDataset object.
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.PixelRepresentation = 0
    ds.HighBit = 11
    ds.BitsStored = 12
    ds.BitsAllocated = 16
    ds.SmallestImagePixelValue = '\\x00\\x00'
    ds.LargestImagePixelValue = '\\xff\\xff'

    ds.Columns = pixel_array.shape[0]
    ds.Rows = pixel_array.shape[1]
    ds.is_little_endian = True

    # set the image direction
    sor = np.array(dct_info['space origin'])
    np_sor = sor.astype(np.float)

    # figure out the pixel spacing
    a = np.array(dct_info['space directions'])
    b = a.astype(np.float)
    da = np.diag(b)
    da_inv = 1. / da
    c = b * (np.eye(3) * da_inv)
    ds.PixelSpacing = [str(x) for x in da[:-1]]
    ds.ImageOrientationPatient = [str(x) for x in c.flatten()[:-3]]

    sthick = da[-1]
    ds.SliceThickness = str(sthick)
    Nslice = dct_info['sizes'][2]

    #     slice_loc = str(np_sor[-1] + (Nslice - slice_idx)*sthick)
    slice_loc = str(np_sor[-1] + (slice_idx) * sthick)

    #     print(slice_loc)
    ds.SliceLocation = slice_loc
    imgpospatient = copy.copy(dct_info['space origin'])
    imgpospatient[2] = slice_loc
    #     print(imgpospatient)
    ds.ImagePositionPatient = imgpospatient

    if pixel_array.dtype != np.uint16:
        pixel_array = pixel_array.astype(np.uint16)
    ds.PixelData = pixel_array.T.tostring()

    ds.save_as(filename)
    return

def Nrrd2Dicom(nrrd_fname, write_dir, filetag):
    np_data, info = nrrd.read(nrrd_fname)
    # make direction for dicoms if doesn't exist
    filetag = (os.path.splitext(os.path.basename(nrrd_fname))[0])
    dicom_dir = '{}/{}_dicom'.format(write_dir, filetag)

    try:
        os.makedirs(dicom_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for ii in range(np_data.shape[2]):
        np_array = np_data[:, :, ii]
        fname = '{}/{:0>4}.dcm'.format(dicom_dir, ii)
        array2Dicom(ii + 1, np_array, info, fname)

