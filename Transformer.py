from __future__ import division
import numpy as np
import math
import os
import errno
import pandas as pd
import re
import shutil
from collections import OrderedDict
import os.path

"""
note: don't need get element composition info as it's not really used in the Data.dat in VHDMSD not everything needs to in a class in python
@staticmethod is mostly useless in python

# 12.18.2014:
    - When incrementing steps to build an array with non-integer step, use np.linspace instead of np.arange for consistent result
    - see this post: http://stackoverflow.com/questions/10011302/python-numpy-arange-unexpected-results
    - nn = np.arange(-halfn*dn,dn*(halfn+1),dn,dtype='float32') output inconsistent output, especially with float32 instead of float64
    - much more consistent when using np.linspace; the below linspace statement in GetDimension(n,dn) has the same output as matlab statement -(dx*halfnx):dx:(dx*halfnx) etc.
    - this may be the reason why the Geant4 VHDMSDv3 kept getting stuck particles between logical volume 'RepX' and 'phantom'
"""

def IsEven(x):
# a helper method; does not need to be in a class: return true if x is even, otherwise
    return x % 2 == 0

def GetDimension(n,dn):
#a helper method; does not need to be in a class
    if IsEven(n):
        halfn = n/2.
        nn = np.linspace(-halfn*dn,halfn*dn,n+1)
    else:
        halfn = math.ceil(n/2.)
        halfdn = dn/2.
        nn = np.linspace(-halfn*dn+halfdn,halfn*dn-halfdn,n+1)
    return nn

def MakeDir(fdir):
# source: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary
    try:
        os.makedirs(fdir)
        #os.mkdir(fdir) # only makes a directory
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print '\nBe careful! directory %s already exists!' % fdir
            
            
def IsInList(a,alist):
    for t in alist:
        if re.search(t,a):
            return True
    return False
    
def ListIntersect(a,b):
    """
    return a list with length of b where true if a == b else false
    for all elements in a, find the first element in b that match a
    """
    output = [False]*len(b)
    for s in a:
        for i in range(len(b)):
            if s == 'Cranium' and s == b[i]:  #exception when matching for 'cranium'
                output[i] = True
                break
            elif s != 'Cranium' and (re.search(re.escape(s),b[i]) or re.search(b[i],s)):  # make sure the pattern in escape character with string such as 'xxx([' containing non-letter characters
                output[i] = True
                break
    return output
    

class Transformer:
    def __init__(self,*args,**kwargs):
        self.theFtype,self.theFrdatadir,self.thePTdir,self.theFwdir,self.geotag,self.ecomptag = args[0:]
        
        self.theFwgeodir = '{}/GeometryIM/G4IM/{}'.format(self.theFwdir,self.geotag)
        self.theFwsrcdir = '{}/SourceMap/{}'.format(self.theFwdir,self.geotag)
        binfwdir = '{}/GeometryIM/binIM/{}'.format(self.theFwdir,self.geotag)
        self.theVolfile = '{}/GeoVol.bin'.format(binfwdir)
        
        #if kwargs.get('flatvol') != None:  # case: no G4 Geometry or SourceMap files exist yet
        if isinstance(kwargs.get('flatvol'), np.ndarray):
            self.flatvol = kwargs.get('flatvol')

            # create a directory for file writing for G4IM and binIM files
            MakeDir(self.theFwgeodir)
            MakeDir(binfwdir)
        else:  # GeoVol.bin exists and be sure to load it in the class when needed
            self.flatvol = None

        if kwargs.get('nxyz'):
            self.nx,self.ny,self.nz = kwargs.get('nxyz')
        else:
            self.nx = self.ny = self.nz = None
        if kwargs.get('dxyz'):
            self.dx,self.dy,self.dz = kwargs.get('dxyz')
            self.dxyz = 0.001*self.dx*self.dy*self.dz  # unit: cm^3
        else:
            self.dx = self.dy = self.dz = None
        
        self.thex0 = self.thex1 = self.they0 = self.they1 = 0
        self.thezz = np.empty(0)
        self.theVol = np.empty(0,dtype=self.theFtype)
        self.theElecomp = np.empty(0)
        self.theOrganInfo = {}
        self.theEcompfname = '{}/OrganData/ECompDensity_ufh{}.txt'.format(self.theFrdatadir,self.ecomptag)
        self.theOrganNamefname = '{}/OrganData/OrgantagvsName_general.txt'.format(self.theFrdatadir)
        
        
        
    def WriteGeoBin(self):
        flatArray = np.ravel(self.flatvol).astype(self.theFtype)
        flatArray.tofile(self.theVolfile)
        print 'write geo bin: {}'.format(self.theVolfile)
        
    def SetTheVol(self):
        # make sure the volume is in the right format
        self.theVol = self.flatvol.reshape((self.nz,self.ny,self.nx),order='C')
    

    def Bin2Vol(self):      #instance member function
        tmp = np.fromfile(self.theVolfile,dtype=self.theFtype)
        self.theVol = tmp.reshape((self.nz,self.ny,self.nx),order='C')

           
    def GetOrganInfo(self):  #instance member function
        # read organtag vs name info from fname into a dictionary
        # copy the fname to the G4 geometry dir
        tmp = {}
        fname = '{}/OrgantagvsName.txt'.format(self.theFwgeodir)
        if not os.path.isfile(fname):
            #print 'no found set path..'
            tmpfname = '{}/GeoIMdata/{}/OrgantagvsName.txt'.format(self.thePTdir,self.geotag)
            # copy the organtag vs name file to geodir
            shutil.copy2(tmpfname,self.theFwgeodir)

        with open(fname,'r') as text:
            for line in text:
                tag,name = line.split()
                if tmp.has_key(tag):
                    print 'found duplicate! organ name: {}'.format(tag)
                else:
                    tmp[name] = int(tag)
                    #print name,tag
                    
       
        
        # assigned the sorted dictionary (by value) to a class variable
        self.theOrganInfo = OrderedDict(sorted(tmp.items(),key=lambda x: x[1]))  # sort by value, return a OrderedDict object, which is a dictionary
        #self.theOrganInfo = sorted(tmp.items(),key=lambda x: x[1])  # sort by value, this return a list
        
        # sort the dictionary and return list of tuples
        #self.theOrganInfo = sorted(tmp.items())  # sort by key
        #self.theOrganInfo = sorted(tmp.items(),key=lambda x: x[1])  # sort by value
    
    def GetECompInfo(self):
        # get the organ info and create the appropriate elemental composition file for simulation
        # for now, the tissue composition is based on 10-yr UFH human phantom
        Ecomp_base = pd.read_csv(self.theEcompfname,sep='\t',skiprows=1,header=None)
        Organname_base = pd.read_csv(self.theOrganNamefname,sep='\t',header=None)
       
        name = self.theOrganInfo.keys()
        tag = self.theOrganInfo.values()
        theEcomp = np.empty([len(name),Ecomp_base.shape[1]],dtype=object)
        thelcOrganname_base = Organname_base[1].str.lower()
        for i in range(len(name)):
            thelcname = name[i].lower()
            if thelcname in ['cranium','bone']:  # match word by word for 'Cranium'
                theEcomp[i,1:] = Ecomp_base.ix[Organname_base[1] == 'Cranium',1:].as_matrix()  #.ix is for mixed integer and lable based access vs. .loc is strickly label based access
                #print 'i = {},organ: {}, matched organ: {}'.format(i,name[i],c)
                #print Ecomp_base.ix[Organname_base[1] == c,1:].as_matrix()[0]
            elif thelcname == 'sp-cranium':
                theEcomp[i,1:] = Ecomp_base.ix[thelcOrganname_base == thelcname,1:].as_matrix()
                #print 'i = {},organ: {}, matched organ: {}'.format(i,name[i],c)
                #print Ecomp_base.ix[Organname_base[1] == c,1:].as_matrix()[0]
            else:
                for c in Organname_base[1].tolist():
                    thelcc = c.lower()
                    #if re.search(re.escape(thelcc),thelcname) or re.search(thelcname,thelcc):
                    if re.search(re.escape(thelcc),thelcname) or re.search(thelcc,thelcname) or re.search(thelcname,re.escape(thelcc)) or re.search(thelcname,thelcc):
                        theEcomp[i,1:] = Ecomp_base.ix[thelcOrganname_base == thelcc,1:].as_matrix()[0]  #take the first name found in the organname base
                        #print 'i = {},organ: {}, matched organ: {}'.format(i,name[i],c)
                        #print Ecomp_base.ix[Organname_base[1] == c,1:].as_matrix()[0]
                        break
            
        # assign the correct organ tag
        theEcomp[:,0] = map(int,tag)
        
        # write ecomp files
        fname = '{}/ECompDensity.txt'.format(self.theFwgeodir)
        thefmt = '\t'.join(['%d'] + ['%2.2f']*(theEcomp.shape[1]-2)+['%2.4f'])
        np.savetxt(fname,theEcomp,fmt=thefmt,header=str(theEcomp.shape[0]),comments='')
        
        
    
    def GetMCInfo(self):
        # copy the two energybin files
        for i in range(1,3):
            efname = '{}/OrganData/Energybin{:d}.txt'.format(self.theFrdatadir,i) 
            shutil.copy2(efname,self.theFwgeodir)
        
        # make the organ tag of interest
        # all organs
        fname1 = '{}/OrgantagOfInterest_all.txt'.format(self.theFwgeodir)
        tmp = np.array(self.theOrganInfo.values())
        np.savetxt(fname1,tmp,fmt='%d',header=str(len(self.theOrganInfo)),comments='')
        #f = open(fname1,'w')
        #tmp = [len(self.theOrganInfo)] + self.theOrganInfo.values()
        #f.write('\n'.join(map(str,tmp)))
        #f.close()
        
        # only skeletal organs
        sktag = np.array([],dtype=int)
        for k,val in self.theOrganInfo.items():
            if re.search('Bone',k) or re.search('sp-',k) or re.search('mc-',k) or re.search('Cranium',k):
                sktag = np.append(sktag,val)
        fname2 = '{}/OrgantagOfInterest_skel.txt'.format(self.theFwgeodir)
        np.savetxt(fname2,sktag,fmt='%d',header=str(len(sktag)),comments='')
        
        #f = open(fname2,'w')
        #tmp = [len(sktag)] + sktag
        #f.write('\n'.join(map(str,tmp)))
        #f.close()
            
    
    def GetPixelDim(self):
        # get pixel spacing
        ll = [(self.nx,self.dx),(self.ny,self.dy),(self.nz,self.dz)]
        xx,yy,self.thezz = [GetDimension(t,dt) for t,dt in ll]
        self.thex0 = xx[0]
        self.thex1 = xx[-1]
        self.they0 = yy[0]
        self.they1 = yy[-1]
        
    def Vol2G4file(self):  #instance member function
        print 'Transformer: start to write .g4m files..'
        # write volume and organtag info to a .g4m
        fnamelist = []
        for i in range(self.nz):
            ftag = 'G4_{:0>4d}.g4m'.format(i)
            fname = '{}/{}'.format(self.theFwgeodir,ftag)
            fnamelist.append(ftag)
            f = open(fname,'w')  # check if fname exists?
            # write organ info
            f.write('%d\n' % len(self.theOrganInfo.items()))
            for k,val in sorted(self.theOrganInfo.items(),key=lambda x: x[1]):  # sort the dictionary by key x[1], sorty by value x[0]
                f.write('%s %s\n' % (val,k))
            
            # write image slice x,y,z info
            f.write('%d %d %d\n%4.3f %4.3f\n%4.3f %4.3f\n%4.3f %4.3f\n' % (self.nx,self.ny,1,self.thex0,self.thex1,self.they0,self.they1,self.thezz[i],self.thezz[i+1]))
            
            # write the pixel value
            imtuple = tuple(map(tuple,self.theVol[i,:,:]))
            f.write('\n'.join(" ".join(map(str, x)) for x in imtuple))  # map function: apply *the function* over items of iterable and return a list of result.  
            #e.g. map(str,x): apply str(n) for all elements in x and join them together with ' '
            f.close()
        
        # write Data.dat
        fname = '%s/Data.dat' % (self.theFwgeodir)
        f = open(fname,'w')
        f.write('%d\n' % len(fnamelist))
        f.write('\n'.join(ii for ii in fnamelist))
        f.close()
        
    def makeSourceMapFile(self,srcname):
        # easy case: a source organ defined by one organ name
        # tricky case: a source organ can be defined by several organ names in theOrganInfo, srcname can be [('xxx','iii'),'yyy',('zzz','www')]
        # solution: input srcname is a dictionary with {'src1':['xx','yy' ...], 'src2':['xx','zz',...]..} for all source organs that may include multiple organ tags
        
        # find the organ tag(s) for a given src organ
        for ss in srcname.keys():
            srclist = [self.theOrganInfo[t] for t in srcname[ss]]
            isSrc = np.in1d(self.theVol,srclist)
            if np.any(isSrc):  # make sure there's some intersection between theVol and srclist, np.all is much FASTER than sum()
                probmap = np.zeros(self.nx*self.ny*self.nz,dtype='float32')
                # a boolean array (same size as theVol) where True when theVol[i] is in srclist False otherwise
                # make sure srclist is a list not a numpy array (huge time suck when it's a numpy array)
                probmap[isSrc] = 1.0
                probmap = probmap/np.sum(probmap)

                # sparsify the array
                indx = np.nonzero(probmap != 0)[0]
                val = probmap[indx]
                output = np.hstack((indx[:,np.newaxis],val[:,np.newaxis]))
                
                # write to files & create a directory for file writing
                thesrcdir = '{}/{}'.format(self.theFwsrcdir,ss)
                MakeDir(thesrcdir)
                fname = '{}/SparseDoseMap.g4d'.format(thesrcdir)
                f = open(fname,'w')  # check if fname exists?
                # write image slice x,y,z info
                f.write('%d %d %d\n%4.3f %4.3f\n%4.3f %4.3f\n%4.3f %4.3f\n' % (self.nx,self.ny,self.nz,self.thex0,self.thex1,self.they0,self.they1,self.thezz[0],self.thezz[-1]))
                outputtuple = tuple(map(tuple,output))
                f.write('\n'.join('{:.0f} {:2.5e}'.format(x[0],x[1]) for x in outputtuple))  # print the (indx,val) tuple in the desired print format
                f.close()
                print 'Transformer: write sourcemap: {}!'.format(fname)
            else:
                print 'Transfomer: the specificed organtag were not found in the geometry volume!'

    def ComputeOrganMass(self,srcname):
        """
            Compute the mass of all segmented organs in the geometry volume
        """
        # load in the elemental composition
        ecompfname = '{}/ECompDensity.txt'.format(self.theFwgeodir)
        ecompdf = pd.read_csv(ecompfname,sep='\t',skiprows=1,header=None)
        
        elenamefname = '{}/OrganData/ElementAtomicN.txt'.format(self.theFrdatadir)
        elenamedf = pd.read_csv(elenamefname,sep='\t',skiprows=0,header=None)
        ecompdf.columns = ['OrganTag'] + elenamedf.iloc[:,0].tolist() + ['Density (g/cm3)']
        
        # compute mass and volume
        self.theOrganMassDF = pd.DataFrame(columns=['OrganName','Volume (cm3)','Mass (g)'],index=range(len(srcname)))
        combo = zip(range(len(srcname)),srcname.items())
        for ii,oo in combo:
            mass = 0.0
            vol = 0.0
            for oname in oo[1]:  # go through all elements in the organ tag(s) that define a given srcname
                otag = self.theOrganInfo[oname]
                isOrgan = np.in1d(self.theVol,otag)
                # rho = ecompdf.loc[ecompdf['OrganTag'] == otag,'Density (g/cm3)'].iget(0,axis=1)
                rho = ecompdf.loc[ecompdf['OrganTag'] == otag, 'Density (g/cm3)'].iat[0]
                if np.any(isOrgan):
                    vol = vol + np.sum(isOrgan)*self.dxyz  #unit: cm^3
                    mass = mass + vol*rho  #unit: g
            tmp = pd.Series({'OrganName':oo[0],'Volume (cm3)':vol,'Mass (g)':mass})
            self.theOrganMassDF.ix[ii] = tmp
        
    """   
    def ComputeOrganMass(self):
        #Compute the mass of all segmented organs in the geometry volume
        
        # load in the elemental composition
        ecompfname = '{}/ECompDensity.txt'.format(self.theFwgeodir)
        ecompdf = pd.read_csv(ecompfname,sep='\t',skiprows=1,header=None)
        
        elenamefname = '{}/OrganData/ElementAtomicN.txt'.format(self.theFrdatadir)
        elenamedf = pd.read_csv(elenamefname,sep='\t',skiprows=0,header=None)
        ecompdf.columns = ['OrganTag'] + elenamedf.iloc[:,0].tolist() + ['Density (g/cm3)']
        
        # compute mass and volume
        self.theOrganMassDF = pd.DataFrame(columns=['OrganName','Volume (cm3)','Mass (g)'],index=range(len(self.theOrganInfo)+1))
        combo = zip(range(len(self.theOrganInfo)),self.theOrganInfo.items())
        for ii,oo in combo:
            isOrgan = np.in1d(self.theVol,oo[1])
            rho = ecompdf.loc[ecompdf['OrganTag'] == oo[1],'Density (g/cm3)'].iget(0,axis=1)

            if np.any(isOrgan):
                vol = np.sum(isOrgan)*self.dxyz  #unit: cm^3
                mass = vol*rho  #unit: g
                tmp = pd.Series({'OrganName':oo[0],'Volume (cm3)':vol,'Mass (g)':mass})
                self.theOrganMassDF.ix[ii] = tmp
        # fill in 'totalbody'
        tmp = pd.Series({'OrganName':'TotalBody','Volume (cm3)':self.theOrganMassDF['Volume (cm3)'].sum(),'Mass (g)':self.theOrganMassDF['Mass (g)'].sum()})
        self.theOrganMassDF.ix[ii+1] = tmp
    """

                
                





