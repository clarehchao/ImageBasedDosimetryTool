#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 3/19/18

@author: shuang
@goal: convert NRRD file from 3dslicer to dicom to import into Amide for further VOI editting

"""

from ImVolTool import FileIO as fi
import sys
import os



if __name__ == '__main__':
    finput = sys.argv[1]
    param_dict = fi.LoadInputParameter(finput)

    # Set up MIBG patient data directory
    ptdir = '{}/PT{:0>4d}'.format(param_dict['NBptdir'],param_dict['PTid'])
    slicer_dir = '{}/3dslicer'.format(ptdir)
    lst_vois = param_dict['organvoiname']

    # write dicoms for all the VOIs for a given subject
    for oo in lst_vois:
        print('voi name: {}'.format(oo))
        nrrd_fname = '{}/{}-label.nrrd'.format(slicer_dir, oo)
        filetag = (os.path.splitext(os.path.basename(nrrd_fname))[0])
        fi.Nrrd2Dicom(nrrd_fname, slicer_dir, filetag)