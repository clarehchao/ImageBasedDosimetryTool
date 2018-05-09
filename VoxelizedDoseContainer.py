from __future__ import division
import numpy as np
from ImVolTool import xyz2vol as xv
import os
import errno
import re

"""
NOTE:
11/21/2014: somehow loading the Edep vol from Coord2Vol from xyz2vol doesn't work...? why???
            Coord2Vol(xxx) doesn't seem to work properly called from this class definition....

"""


def MakeDir(fdir):
# source: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary
    try:
        os.mkdir(fdir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print '\nBe careful! directory %s already exists!' % fdir

class VoxelizedDoseContainer:
    def __init__(self,*args):
        # arguments for the class constructor
        self.theGeoDir = args[0]
        self.theGeoTag = args[1]
        self.theNxyz = args[2]
        self.theDxyz = args[3]
        self.theNxNyNz = np.prod(np.array(self.theNxyz))
        
        # class private variables
        self.theDoseDir = ''
        self.theGeoVol = np.empty(0,dtype='uint8')
        self.theinfodict = {}
        self.edepftype = np.float32
        
        # add more dose volume e.g. fluence volumn when needed
        self.constK = 2.13*10/(3.7e-2)/3600 # unit: mGy-g/(MBq-s-MeV)
        self.voxvol = np.prod(0.1*np.array(self.theDxyz))  # voxel volume, unit: cm^3
        #self.Gy2MeVperGram = (1/(1.60217646e-16))*1.0e-3/1000
        #organname = ['Bone','Lungs','Adrenal','Brain','Heart','Liver','SalivaryGlands','Spleen']
        #tumorname = ['Tumor%d' % n for n in range(1,13)]
        #self.TargetOrgan = organname + tumorname
    
    def getInfo2Dict(self):
        # construct a dictionary of organ name storing organ tag and material density
        fname1 = '{}/G4IM/{}/ECompDensity.txt'.format(self.theGeoDir,self.theGeoTag)
        data1 = np.loadtxt(fname1,skiprows=1)
        tmpdict = {}
        for i in range(data1.shape[0]):
            tmpdict[int(data1[i,0])] = data1[i,-1]
        
        # read a text file and save to a dictionary
        fname2 = '{}/G4IM/{}/OrgantagvsName.txt'.format(self.theGeoDir,self.theGeoTag)
        f = open(fname2,'r')
        for line in f:
            splitline = line.split()
            organtag = int(splitline[0])
            name = splitline[1]
            if name not in self.theinfodict:
                self.theinfodict.setdefault(name,[])
            self.theinfodict[name].append(organtag)
            self.theinfodict[name].append(tmpdict[organtag])
        
        # define target organ
        self.TargetOrgan = [s for s in self.theinfodict.keys() if not (re.search('Air',s) or re.search('Teeth',s))]
    
    def loadGeoVol(self):
        thefile = '{}/binIM/{}/GeoVol.bin'.format(self.theGeoDir,self.theGeoTag)
        print thefile,self.theGeoDir
        self.theGeoVol = xv.FlattenFile2Vol(thefile,'uint8',self.theNxyz[0],self.theNxyz[1],self.theNxyz[2])
    
    def loadEdepVol(self,fdir,run1,run2):
        # figure out how to pass in  $d in string with run1 and run2 for dose accumulation!
        self.theDoseDir = fdir
        self.simNevent = 0  #update the number of events in the simulation for dose output
        self.theEdepVol = np.zeros(self.theNxNyNz,dtype=np.float32)
        for r in range(run1,run2):
            #thefile = '{}/Run{}/BustOutRoot/Edep.dat'.format(self.theDoseDir,r)
            thefile = '{}/Run{}/RootData/Edep.root'.format(self.theDoseDir, r)
            tmp = xv.Coord2Vol(thefile,self.edepftype,self.theNxyz)
            self.theEdepVol = self.theEdepVol + tmp
        
            # Determine the number of total events in the dose result file
            thelogfile = '{}/Run{}/log.txt'.format(self.theDoseDir,run1)
            f = open(thelogfile,'r')
            self.simNevent = self.simNevent + int(re.findall(r'Number of Events in this run: ([\w]+)',f.read())[0])
            
    def loadOtherEdepVol(self,fdir,ftag,run1,run2):
        # figure out how to pass in  $d in string with run1 and run2 for dose accumulation!
        self.theDoseDir = fdir
        self.theOtherVol = np.zeros(self.theNxNyNz,dtype=np.float32)
        for r in range(run1,run2):
            thefile = '{}/Run{}/BustOutRoot/{}Edep.dat'.format(self.theDoseDir,r,ftag)
            tmp = xv.Coord2Vol(thefile,self.edepftype,self.theNxyz)
            self.theOtherVol = self.theOtherVol + tmp
        
    def ComputeSvalue(self):
        """
        :return: S-value is in unit of mGy/(MBq-s)
        """
        #print 'In ComputeSvalue():',self.theinfodict
        
        # empty out the dictionary before the next src organ
        self.theSVdict = {}
        for oo in self.TargetOrgan:
            tag = self.theinfodict[oo][0]
            indx = (self.theGeoVol == tag)
            totEdep = np.sum(self.theEdepVol[indx])
            mass = self.voxvol*self.theinfodict[oo][1]*np.sum(indx)
            self.theSVdict[oo] = self.constK*totEdep/mass/self.simNevent
    
    def ComputeSvalueStats(self,dosedir,runlist,drun):
        # Energy deposit compilation
        edeptmp = np.zeros((len(self.TargetOrgan)+1,len(runlist)),dtype=np.float32)
        for ii in range(len(runlist)):
            self.loadEdepVol(dosedir,runlist[ii],runlist[ii]+drun)
            self.ComputeSvalue()
            edeptmp[:,ii] = np.array([self.simNevent] + self.theSVdict.values())

        # construct final Edep w/ stats dictionary
        edepMean = np.mean(edeptmp,axis=1)
        edepStd = np.std(edeptmp,axis=1)
        organname = self.theSVdict.keys()
        
        # empty out the dictionary before the next src organ
        self.edepstats = {}
        for ii in range(len(organname)):
            k = organname[ii]
            if k not in self.edepstats:
                self.edepstats.setdefault(k,[])
            self.edepstats[k].append(edepMean[ii+1])
            self.edepstats[k].append(edepStd[ii+1])
        self.SimfinalNevent = edepMean[0]


# need a function to pull out from mysql DB to look at the data of interest
# e.g. compare the dose result (a src-to-target organs) from different simulation pkgs
# pull out mass, position related information of organs for different src-to-target organ pairs
# other mysql request to get data out of DB to present for interesting visualizations
                
            
            
            
            
            
    
    
