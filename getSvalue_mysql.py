#! /usr/bin/env python
# -*- coding: utf-8 -*-

from DBTool import DoseDB as ddb
from ImVolTool import xyz2vol as xv
from ImVolTool import FileIO as fi
import numpy as np
import VoxelizedDoseContainer as vdc
import MySQLdb as mdb
import os
import sys

if __name__ == '__main__':
    finput = sys.argv[1]
    param_dict = fi.LoadInputParameter(finput)
    
    # set up data directories
    geodir = '{}/G4INPUT/GeometryIM'.format(param_dict['G4dir'],param_dict['geotag'])
    datadir = '{}/G4.9.6.p02work/{}/data/GEO_{}/SRCMP_{}'.format(param_dict['G4dir'],param_dict['G4AppName'],param_dict['geotag'],param_dict['geotag'])
    
    # find all src organ directories
    #lsubdir = os.listdir(datadir)
    #srconame = [s for s in lsubdir if os.path.isdir('{}/{}'.format(datadir,s)) and s not in param_dict['excludesrcname']]
    srconame = param_dict['srcname'].keys()

    # initialize an instance of vdc
    thedoseobj = vdc.VoxelizedDoseContainer(geodir,param_dict['geotag'],param_dict['nxyz'],param_dict['dxyz'])
    
    # Set up organ and geometry info for the vdc instance
    thedoseobj.getInfo2Dict()
    thedoseobj.loadGeoVol()
    
    # start the connection with database and others
    # this works on Higgs as if one specifies 'localhost', it tries to connect via default /tmp/mysql.socket
    # see http://stackoverflow.com/questions/4662364/cant-connect-to-localhost-using-pythons-mysqldb
    #con = mdb.connect('localhost','testuser','test000','UCSFDoseDB')
    con = mdb.connect('127.0.0.1','testuser','test000','UCSFDoseDB')
    simpkg = 'G4.9.6.p02'

    # check SimInfo and DoseInfo tables exist
    if ddb.CheckTableExist(con,'SimInfo') and ddb.CheckTableExist(con,'DoseInfo'):
        print 'SimInfo and DoseInfo tables exist in the database!'
    else:
        'SimInfo and DoseInfo tables do not exist in the database! Create the tables!'
        ddb.CreateTableDB_DoseSimInfo(con)
    
    # go through all the source organs for a given geo setup
    for ss in srconame:
        print 'source organ: {}'.format(ss)
        dosedir = '{}/{}/{}'.format(datadir,ss,param_dict['srcparticle'])
       
        # figure out how to bundle the MC Run's for Statistics
        #lrundir = [s for s in os.listdir(dosedir) if os.path.isdir('{}/{}'.format(dosedir,s))]
        #run1,run2 = [1,len(lrundir)]
        #runlist = range(run1,run2+1,param_dict['drun'])
        runlist = range(param_dict['therun'][0],param_dict['therun'][1]+1,param_dict['therun'][2])

        # Energy deposit statistics compilation
        thedoseobj.ComputeSvalueStats(dosedir,runlist,param_dict['therun'][2])
        statsdict = thedoseobj.edepstats
        
        # Print out the Svalue stats result
        for k,val in statsdict.items():
            print 'k = {}, val = {}'.format(k,val)
        print 'Simulation final Nevent: {}'.format(thedoseobj.SimfinalNevent)

        # Update database with the new calculation
        varags = [simpkg,param_dict['geotag'],param_dict['srcparticle'],ss,thedoseobj.SimfinalNevent,len(runlist),param_dict['G4AppName'],param_dict['G4SimDataDir']]
        ddb.Insert2DB_fancydict_DoseSimInfo(con,statsdict,varags)







