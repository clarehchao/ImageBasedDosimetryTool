#! /usr/bin/env python
# -*- coding: utf-8 -*-

#import the toolbox for converting a flat coor txt file to a flat or 3d volum
from ImVolTool import xyz2vol as xv
from ImVolTool import FileIO as fi
import Transformer as tf
import sys
from DBTool import DoseDB as ddb
import MySQLdb as mdb

"""
   Objective:
   - Load in the necessary input parameters for a given geometry for G4 MC
   - Call an instance of Transformer class given the GeoVol.bin already exists
   - Write additional source map files without re-writing G4IM and binIM files again
"""

if __name__ == '__main__':
    finput = sys.argv[1]
    param_dict = fi.LoadInputParameter(finput)
    
    # Set up MIBG patient data directory
    ptdir = '{}/PT{:0>4d}'.format(param_dict['NBptdir'],param_dict['PTid'])

    # convert the segmented volume to G4-friendly files
    tf = tf.Transformer(param_dict['binfwtype'],param_dict['VHDMSDdir'],ptdir,param_dict['fwdir'],param_dict['geotag'],param_dict['ecomptag'],nxyz=param_dict['nxyz'],dxyz=param_dict['dxyz'])

    # convert a binary file to a volume
    tf.Bin2Vol()
    
    # get the appropriate info about the volume for G4file writing
    tf.GetPixelDim()
    tf.GetOrganInfo()
    
    tf.ComputeOrganMass(param_dict['srcname'])
    print tf.theOrganMassDF
    

    # Insert the mass info into the mysql database
    db_pw = ddb.get_DB_auth_info(param_dict['DB_auth_dir'], param_dict['DB_usr'])
    con = mdb.connect(host='127.0.0.1', user=param_dict['DB_usr'], passwd=db_pw, db='UCSFDoseDB')

    # check to see if the GeoInfo table exists
    tablename = 'GeoInfo'
    if ddb.CheckTableExist(con,tablename):
        print 'table exist! no need to create one!'
    else:
        print 'table does not exist!'
        ddb.CreateTableDB_GeoInfo(con)
    varags = param_dict['geotag']
    ddb.Insert2DB_GeoInfo(con,tf.theOrganMassDF,varags)
