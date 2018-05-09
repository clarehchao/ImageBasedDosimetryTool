from __future__ import division
import numpy as np
import pydicom  # dicom reader package
import glob   # get all the files in a directory
import sys
import os
import re
import errno
import pandas as pd
import ROOT
from root_numpy import tree2array

def IsInBound(val,low,high):
    """
    check to make sure VAL in within the range (low,high]
    """
    if len(val) == 1:   # for a single value
        if val < low:
            #print 'less than low'
            return low
        elif val >= high:
            #print 'more than high'
            return high-1
        else:
            return val
    elif len(val) > 1:  # for the case of an array
        val[val < low] = low
        val[val >= high] = high-1
        return val
    else:
        print 'empty input!'

def properRound(x):
    """
    Numpy rounds x.5 to th nearest even number (not sure why...)
    this function does the proper round e.g. for x.5, round up to the nearest integer
    Input: x is a single scalar or a numpy array
    Output: return a scalar or an array of properly rounded numbers
    """
    y = x - np.floor(x)
    if len(x) == 1:  # a scalar
        if (0 < y < 0.5):
            return np.floor(x)
        else:
            return np.ceil(x)
    else:  # an array
        output = np.empty(x.shape,dtype='int')  # make sure the output is in 'int' not 'float'
        indx = np.all([y > 0,y < 0.5],axis=0)
        notindx = np.logical_not(indx)
        output[indx] = np.floor(x[indx],out=np.empty_like(x[indx],dtype=np.int_),casting='unsafe')
        output[notindx] = np.ceil(x[notindx],out=np.empty_like(x[notindx],dtype=np.int_),casting='unsafe')
        return output
        
def GetUniqueCount(x):
    # input: x is a numpy array
    thecount = {}
    for v in x:
        if v not in thecount:
            thecount[v] = 1
        else:
            thecount[v] += 1
    return thecount
    
def SaveFlattenVol(vol,fname):
    outputType = 'float32'
    outputExtension = '.f32le'
    lefile = fname + outputExtension
    #flatArray = vol.flatten().astype(outputType)
    flatArray = np.ravel(vol).astype(outputType)
    flatArray.tofile(lefile)
    
    
def Xyz2Vol(fname,ftype,dxyz,xyz0,nxyz,nskiprow,dataindx):
    # get the data into workable form
    #data = np.loadtxt(fname,skiprows=1)
    #val,wt,xm,ym,zm = np.hsplit(data,data.shape[1])
    #ix0,iy0,iz0,xm,ym,zm = np.hsplit(data,data.shape[1])

    # a more general approach to parse the imported data
    data = np.loadtxt(fname,skiprows=nskiprow)
    xm,ym,zm = np.hsplit(data[:,dataindx],len(dataindx))
    #print np.amin(xm),np.amax(xm)
    #print np.amin(ym),np.amax(ym)
    #print np.amin(zm),np.amax(zm)

    # get dxyz and xyz dimension
    dx,dy,dz = dxyz
    x0,y0,z0 = xyz0
    #print nxyz
    nx,ny,nz = nxyz
    #ixx = IsInBound(properRound((xm-x0)/dx),0,nx)
    #iyy = IsInBound(properRound((ym-y0)/dy),0,ny)
    #izz = IsInBound(properRound((zm-z0)/dz),0,nz)
    
    #tmp = (xm-x0)/dx
    #print tmp
    #print np.round(tmp).astype(np.int)
    ixx = IsInBound(np.round((xm-x0)/dx).astype(np.int),0,nx)
    iyy = IsInBound(np.round((ym-y0)/dy).astype(np.int),0,ny)
    izz = IsInBound(np.round((zm-z0)/dz).astype(np.int),0,nz)
    
    #get axis space
    #xx = np.arange(x0,x0+nx*dx,dx)
    #yy = np.arange(y0,y0+ny*dy,dy)
    #zz = np.arange(z0,z0+nz*dz,dz)
    
    # convert coordinates into indices for 1D array operation
    xyzcoord = np.vstack((izz.T,iyy.T,ixx.T))
    indx = np.ravel_multi_index(xyzcoord,(nz,ny,nx))
    voxct = GetUniqueCount(indx)
    vol_flat = np.zeros(nx*ny*nz,dtype=ftype)
    vol_flat[voxct.keys()] = voxct.values()
    #vol_3d = vol_flat.reshape((nz,ny,nx))
    return vol_flat
    
def Coord2Vol(fname,ftype,nxyz):
    """

    :param fname: the file name of the coordinate text file
    :param ftype: data type for the binary volume
    :param nxyz: dimension of the volume
    :return: a binary volume

    the indexing method produces the same result as Xyz2Vol if one look at vol.reshape((nz,nx,ny))
    """

    if fname.endswith('.dat'):
        data = np.loadtxt(fname)
        print('In xyz2vol:Coord2Vol function, read the .dat file: {}'.format(fname))
    elif fname.endswith('.root'):
        rfile = ROOT.TFile(fname)
        intree = rfile.Get('EdepTree')
        arr = tree2array(intree)
        data = np.stack((arr['posX'], arr['posY'], arr['posZ'], arr['eng']), axis=1)
        print('In xyz2vol:Coord2Vol function, read the .root file: {}'.format(fname))
    else:
        print('::O_O In xyz2vol:Not correct dose data file extension::')

    nx,ny,nz = nxyz
    
    #get the x,y,z position, first three cols of 'data'
    tmp = data[:,:-1]
    #xyz = tmp.conjugate().transpose()
    xyz = tmp.T
    xyzint = xyz.astype(int)
    val = data[:,3]
    
    #vol = np.zeros((nx,ny,nz))
    flatvol = np.zeros(nx*ny*nz,dtype=ftype)
    #convert (x,y,z) indx to flat indices
    indx = np.ravel_multi_index(xyzint,(nx,ny,nz),order='F') # column-major
    #print len(indx)
    #vol.ravel()[indx] = val
    flatvol[indx] = val
    #print np.sum(flatvol)
    return flatvol
    
    
def SaveFlattenVol(vol,fname,ftype):
    #flatArray = vol.flatten().astype(outputType)
    flatArray = np.ravel(vol).astype(ftype)
    flatArray.tofile(fname)
    print 'saved the volume to {} as {}!'.format(fname,ftype)


def VolMask_Threshold(vol,thresh=0):
    # create a vol mask based on a given threshold
    vsize = np.prod(np.array(vol.shape))
    mask = np.zeros(vsize,dtype='uint8')
    if len(thresh) == 2:
        t1,t2 = thresh
        indx = np.all([(vol >= t1),(vol < t2)],axis=0)
        mask[indx] = 1
    else:
        mask[vol.flatten() >= thresh] = 1
    return mask

def SegmentVol_ImMask(vol,maskfname,ftype,segval0=0):
    # create a segmented volume given a list of masks
    vsize = np.prod(np.array(vol.shape))
    #segvol = np.zeros(vsize,dtype=ftype)
    if len(vol.shape) > 1:
        segvol = np.ravel(vol)
    else:
        segvol = vol
    for i in range(len(maskfname)):
        ma = np.fromfile(maskfname[i],dtype=ftype)
        if re.search('Tumor',maskfname[i]) or re.search('lesion',maskfname[i]) or re.search('Kidney',maskfname[i]) or re.search('Heart',maskfname[i]):  # tumor or lesion
            segvol[ma >= 1] = segval0+i
            print 'It\'s a tumor! can overwrite bone voxel!'
        elif re.search('Brain',maskfname[i]): # do not want the image mask to overwrite any 'bone' or 'air' voxel
            indx = np.all([(ma >= 1),(vol != 2),(vol != 0)],axis=0)
            segvol[indx] = segval0+i
            print 'It\'s a Brain, avoid overwriting bone and air voxels!'
            print 'anything above 1? {:f}'.format(np.sum(indx))
        else:
            indx = np.all([(ma >= 1),(vol != 2)],axis=0)
            segvol[indx] = segval0+i
            print 'It\'s Not a brain or tumor or kidney, avoid overwriting bone voxels!'
            print 'anything above 1? {:f}'.format(np.sum(indx))
        
        print 'set a mask value: {},{:f}'.format(maskfname[i],segval0+i)
    return segvol
    
def SegmentVol_ThreshMask(vol,thresh):
    # create a segmented volume given a list of masks
    #vsize = np.prod(np.array(vol.shape))
    segvol = np.zeros(vol.shape,dtype='uint8')
    for i in range(len(thresh)-1):  # go throughe each threshold bin and assign an index of material
        t1 = thresh[i]
        t2 = thresh[i+1]
        indx = np.all([(vol >= t1),(vol < t2)],axis=0)
        #print indx.shape
        segvol[indx] = i
    #print segvol.shape
    return np.ravel(segvol),np.amax(segvol)
    
def FlattenFile2Vol(fname,ftype,nx,ny,nz):
    vol = np.fromfile(fname,dtype=ftype)
    #tmp = np.fromfile(fname,dtype=ftype)
    #vol = tmp.reshape((nz,ny,nx),order='C')
    #vol = tmp.reshape((nz,ny,nx),order='F')
    return vol
    
def Dicom2Vol(dcdir):    #instance member function
    alldcfiles = glob.glob('%s/*.DCM' % dcdir)
    #print alldcfiles[0]
    # find the dcfile tag: e.g. ____I10.DCM
    match = re.search(r'%s/([\w.-]+)I([\w.-]+).DCM' % dcdir,alldcfiles[0])
    if match:
        dcftag = match.group(1)
        #print dcftag
    else:
        raise RuntimeError('Error: cannot get dicom file tag!')
    flist = range(1,len(alldcfiles)+1)
    for i in range(len(flist)):
        dcfile = '%s/%sI%d.DCM' % (dcdir,dcftag,flist[i])
        dc = pydicom.dcmread(dcfile)
        if i == 0:  # initialize
            dxyz = np.array([dc.PixelSpacing[0],dc.PixelSpacing[1],dc.SliceThickness])
            nx = dc.Columns
            ny = dc.Rows
            nz = len(alldcfiles)
            nxyz = np.array([nx,ny,nz])
            vol = np.empty((nz,ny,nx)) 
            #print dc
        im = np.flipud(dc.pixel_array)*dc.RescaleSlope + dc.RescaleIntercept
        #print np.amin(im),np.amax(im)
        vol[i,:,:] = im
    return dxyz,nxyz,vol

def Dicom2Vol_v2(dcdir):    #instance member function
    alldcfiles = sorted(glob.glob('%s/*.dcm' % dcdir))
    
    if alldcfiles:
        for i in range(len(alldcfiles)):
            dcfile = alldcfiles[i]
            dc = pydicom.dcmread(dcfile)
            if i == 0:  # initialize
                dxyz = np.array([dc.PixelSpacing[0],dc.PixelSpacing[1],dc.SliceThickness])
                nx = dc.Columns
                ny = dc.Rows
                nz = len(alldcfiles)
                nxyz = np.array([nx,ny,nz])
                vol = np.empty((nz,ny,nx)) 
                #print dc
            im = np.flipud(dc.pixel_array)*dc.RescaleSlope + dc.RescaleIntercept
            #print np.amin(im),np.amax(im)
            vol[i,:,:] = im
    else:
        print('dicom files are not found!')
            
    return dxyz,nxyz,vol

def MakeDir(fdir):
# source: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary
    try:
        os.mkdir(fdir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print '\nBe careful! directory %s already exists!' % fdir
    

def CombineSegVol(frdir,frtype,baseftag,maskftag):
    # combine the base volume with mask volume(s) and save the volume to the given dir
    # read the base volume
    basefname = '{}/{}.bin'.format(frdir,baseftag)
    vol_base = np.fromfile(basefname,dtype=frtype)
    
    # generate a list of organ mask: in this case, just tumor mask for moby phantoms
    if len(maskftag) > 0:   
        maskfname = ['{}/{}.bin'.format(frdir,n) for n in maskftag]
        maxvval = np.amax(vol_base)
        vol_mask = SegmentVol_ImMask(vol_base,maskfname,frtype,segval0=maxvval+1)
        ovol = vol_mask
        #SaveFlattenVol(vol_mask,gvolfname,fwtype)
    else:
        ovol = vol_base
        #SaveFlattenVol(vol_base,gvolfname,fwtype)
    
    return ovol   


def PMODpix2Vol(fname,nxyz,ftype):
    """
    :param fname: a PMOD .pixeldump filename
    :param nxyz: a numpy array of [nx,ny,nz] wrt the pixel dimension in fname
    :return: a binary mask volume where pixval = 1 if pixel is in fname, 0 otherwise
    """
    df = pd.read_csv(fname,sep='\t',skiprows=1)
    xyzcoord = df.loc[1:,['z','y','x']].as_matrix().astype('int32')
    nx,ny,nz = nxyz
    indx = np.ravel_multi_index(xyzcoord.T,(nz,ny,nx)) # row-major as default order='C'
    #voxct = GetUniqueCount(indx)

    vol_flat = np.zeros(np.prod(nxyz),dtype=ftype)
    #vol_flat[voxct.keys()] = voxct.values()
    vol_flat[np.unique(indx)] = 100
    vol = np.fliplr(np.flipud(vol_flat.reshape((nz,ny,nx))))

    return vol


