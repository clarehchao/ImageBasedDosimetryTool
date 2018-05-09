#! /usr/bin/env python
# -*- coding: utf-8 -*-

#import the toolbox for converting a flat coor txt file to a flat or 3d volum
from ImVolTool import xyz2vol as xv
from ImVolTool import FileIO as fi
import numpy as np
import Transformer as tf
import sys
import os


"""
    Objective:
    if GeoVol.bin does NOT exist, do the following:
    - Load in the necessary input parameters for a given geometry for G4 MC
    - Segment a set of CT images into the organs of interest based on HU thresholds or manual-drawn VOIs from Amide
    - Call an instance of Transformer class given all G4IM and binIM files DO NOT EXIST
    - Write the source map files for the source organs of interest
    Otherwise:
    - load the geometry volume (GeoVol.bin)
    - set up the appropriate image dimension and other info
    - generate the needed source map files

"""

if __name__ == '__main__':
    finput = sys.argv[1]
    param_dict = fi.LoadInputParameter(finput)

    # Set up MIBG patient data directory
    ptdir = '{}/PT{:0>4d}'.format(param_dict['NBptdir'],param_dict['PTid'])

    # check to see if xxxx/binIM/xxxx/GeoVol.bin exists
    geobinfname = '{}/GeometryIM/binIM/{}/GeoVol.bin'.format(param_dict['fwdir'],param_dict['geotag'])
    if os.path.isfile(geobinfname):
        print 'WooHoo! GeoVol.bin exist: load the volume and write the source maps!'
        # convert the segmented volume to G4-friendly files
        tf = tf.Transformer(param_dict['binfwtype'],param_dict['VHDMSDdir'],ptdir,param_dict['fwdir'],param_dict['geotag'],param_dict['ecomptag'],nxyz=param_dict['nxyz'],dxyz=param_dict['dxyz'])

        # convert a binary file to a volume
        tf.Bin2Vol()

        # get the appropriate info about the volume for G4file writing
        tf.GetPixelDim()
        tf.GetOrganInfo()

        # write sourcemap based on the geometry volume
        tf.makeSourceMapFile(param_dict['srcname'])
    else:  # build the volume from the raw CT images to GeoVol.bin
        print 'wait up! GeoVol.bin is not found: let\'s start from the groundup to get the geometry volume ready and write the source maps!'
        # basic CT segmentation based on amide's VOIs
        # 'voifname' include all the names of all amide VOIs (include organs and tumors)
        thedir = '{}/GeoIMdata/{}'.format(ptdir,param_dict['geotag'])
        for nn in param_dict['organvoiname']:
            fname = '{0}/roi_raw_data_{{{1}}}.tsv'.format(thedir,nn)
            vol = xv.Xyz2Vol(fname,param_dict['frtype'],param_dict['dxyz'],param_dict['xyz0'],param_dict['nxyz'],0,param_dict['frindx'])
            binfname = '{}/{}.bin'.format(thedir,nn)
            xv.SaveFlattenVol(vol,binfname,param_dict['binfwtype'])

        # get the CT dicom volume
        if any(ff.endswith('.DCM') for ff in os.listdir(param_dict['ctdir'])):
            ctdxyz, ctnxyz, ctvol = xv.Dicom2Vol(param_dict['ctdir'])
        elif any(ff.endswith('.dcm') for ff in os.listdir(param_dict['ctdir'])):
            ctdxyz,ctnxyz,ctvol = xv.Dicom2Vol_v2(param_dict['ctdir'])
        else:
            sys.exit('WAIT! Did not find any .dcm or .DCM files in {}! Exit program ... '.format(param_dict['ctdir']))
        xv.SaveFlattenVol(ctvol, '{}/ctvol.bin'.format(thedir), 'float32')

        # segment via threshold
        segvol1,maxsegval = xv.SegmentVol_ThreshMask(ctvol,param_dict['HUthresh'])
        xv.SaveFlattenVol(segvol1, '{}/segvol1.bin'.format(thedir), 'uint8')

        # generate a list of organ mask
        voiall = param_dict['organvoiname']
        voifname = ['{}/{}.bin'.format(thedir,ss) for ss in voiall]
        segvol2 = xv.SegmentVol_ImMask(segvol1,voifname,param_dict['binfwtype'],segval0=maxsegval+1)
        xv.SaveFlattenVol(segvol2, '{}/segvol2.bin'.format(thedir), 'uint8')

        # prepare the Organtag vs Name files
        thetag = sorted(np.unique(segvol2))  # check this statement!
        thename = param_dict['HUthresh_name'] + voiall
        fname = '{}/OrgantagvsName.txt'.format(thedir)
        f = open(fname,'w')
        for i in range(len(thetag)):  # sort the dictionary by key x[1], sorty by value x[0]
            f.write('%s %s\n' % (thetag[i],thename[i]))
        f.close()
        
        # convert the segmented volume to G4-friendly files
        tf = tf.Transformer(param_dict['binfwtype'],param_dict['VHDMSDdir'],ptdir,param_dict['fwdir'],param_dict['geotag'],param_dict['ecomptag'],flatvol=segvol2,nxyz=param_dict['nxyz'],dxyz=param_dict['dxyz'])

        # write the binary geo file
        tf.WriteGeoBin()


        # For the instance of transformer, set the main volume to segvol2
        tf.SetTheVol()

        # get the appropriate info about the volume for G4file writing
        tf.GetPixelDim()
        tf.GetOrganInfo()

        # write appropriate files for elemental composition for each tissue material and Monte-Carlo related files
        tf.GetECompInfo()
        tf.GetMCInfo()

        # write G4files of the transformed volume: this step takes longer, can avoid this step if G4files are written already
        tf.Vol2G4file()

        # write sourcemap based on the geometry volume
        tf.makeSourceMapFile(param_dict['srcname'])





