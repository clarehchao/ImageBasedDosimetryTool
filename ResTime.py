import numpy as np
import re
import pandas as pd
import math as mt
from ImVolTool import ImVolFit as ivfit
import MySQLdb as mdb
import glob
import os
import pydicom
from datetime import datetime
from pydicom.tag import Tag
from DBTool import DoseDB as ddb

def check_date_format(ss):
    for fmt in ('%Y%m%d%H%M%S.%f', '%Y%m%d%H%M%S'):
        try:
            return datetime.strptime(ss, fmt)
        except ValueError:
            pass
    raise ValueError('no validate date formate found!')

# make the src_organ name and the ResTimeHrDF organ name matches
def Func_Match_Name1(x,reflst):
    # print('x is {}'.format(x))
    x_strip = ''.join(x.split()).lower()
    refdict = dict(zip(reflst,[''.join(s.split()).lower() for s in reflst]))
    for k,val in refdict.items():
        if re.findall(r'tumor',x_strip):
            # if x_strip == 'tumor6':
            #     print('it is tumor6!')
            #     print(re.findall(r'{}\b'.format(x_strip),val))
            #     print(re.findall(r'{}\b'.format(val),x_strip))
            if re.findall(r'{}\b'.format(x_strip),val) or re.findall(r'{}\b'.format(val),x_strip):
                return k
                break
        else:
            rep = {'(':'\(',')':'\)'}
            rep = dict((re.escape(k), v) for k, v in rep.iteritems())
            pattern = re.compile('|'.join(rep.keys()))
            x_strip_mod = pattern.sub(lambda x: rep[re.escape(x.group(0))], x_strip)
            val_mod = pattern.sub(lambda x: rep[re.escape(x.group(0))], val)
            # if x_strip == 'adrenal':
            #     print('it is adrenal!')
            #     print(x_strip_mod, val)
            #     print(re.findall(r'{}'.format(x_strip_mod),val))
            #     print(val_mod, x_strip)
            #     print(re.findall(r'{}'.format(val_mod),x_strip))
            if re.findall(r'{}'.format(x_strip_mod),val) or re.findall(r'{}'.format(val_mod),x_strip):
                return k
                break

def PETSUV_func(row):
    # print(row['OrganName'])
    if not np.isnan(row['Mean']):
        # print('has mean!')
        return row['Mean']*row['PETSUV_factor']
    elif not np.isnan(row['Averaged[kBq/cc]']):
        # print('has averaged...')
        return row['Averaged[kBq/cc]']*1000.*row['PETSUV_factor']
    else:
        return np.nan

def upperfirst(x):
    return x[:1].upper() + x[1:]

def Func_Match_Name2(x, reflst): # cross match the name of string x with all the strings in reflst
    x = x.lower()
    for s in reflst:
        tmp = s.replace(' ','')
        if x.find(tmp) != -1 or tmp.find(x) != -1:
            # print('x={}'.format(x))
            # print('tmp={}'.format(tmp))
            return s
    # print('no match!')
    return x


class ResTime(object):
    def __init__(self,*args,**kwargs):
        # arguments for the class constructor
        self.PTdir = args[0]
        self.PETsubdirs = args[1]
        self.therapy_isotope_halflife = args[2]  # unit: hour
        self.img_isotope_halflife = args[3]     # unit: hour
        self.srcname_dict = args[4]
        self.Frdatadir = args[5]
        self.DB_auth_dir = args[6]
        self.DB_usr = args[7]

        self.PMODdir = None
        self.theTumordf = None
        self.theOrganDosemGy = None
        self.theResTimeHrDF = None
        
        
        # class private variables
        self.theRTdf = pd.DataFrame()
        
        # therapy Isotope lambda_physical, unit: 1/hour
        self.therapy_isotope_lambda_p = mt.log(2.)/self.therapy_isotope_halflife #unit: 1/hour

        # Imaging Isotope lambda_physical, unit: 1/hour
        self.img_isotope_lambda_p = mt.log(2.)/self.img_isotope_halflife
    
        # the final residence estimate for a given organ and total body with excretion model
        #self.theResTimeHrDF = pd.DataFrame(columns = ['OrganName','p1','p2','p3','p4','r2','Residence Time (Bq-hr/Bq)'])
    
        # mysql connection setup to connect to the local mysql database
        db_pw = ddb.get_DB_auth_info(self.DB_auth_dir, self.DB_usr)
        db_name = 'UCSFDoseDB'
        self.mysqlcon = mdb.connect(host='127.0.0.1', user=self.DB_usr, passwd=db_pw, db=db_name)

        self.TP_slope_idx = [(1,2),(2,3),(1,3)]

    def InitialSetUp(self):
        # Read the dicom header of the IM images to determine time since injection
        self.PETdtSinceInj0 = []
        self.PETSUV_factor0 = []
        for ii in range(len(self.PETsubdirs)):
            petdir = '{}/IM/{}'.format(self.PTdir,self.PETsubdirs[ii])
            # print('petdir = {}'.format(petdir))
            allfiles = os.walk(petdir).next()[2]
            thedc = pydicom.dcmread('{}/{}'.format(petdir,allfiles[0]))

            if thedc.has_key(Tag(0x0009,0x100d)):  # Discovery STE format
                scanDT = check_date_format(thedc[0x0009, 0x100d].value)
            elif thedc.has_key(Tag(0x0008, 0x002a)):  #Philip PET scanner format
                scanDT = check_date_format(thedc[0x0008, 0x002a].value)
            else:
                print("ERROR: scan date dicom tag is INCORRECT!")
           
            tag_Radiopharmaceutical_Information_Sequence = Tag(0x0054, 0x0016)
            tag_Radiopharmaceutical_start_datetime = Tag(0x0018, 0x1078)
            rad_pharm_seq = thedc[tag_Radiopharmaceutical_Information_Sequence]
            if rad_pharm_seq[0].has_key(tag_Radiopharmaceutical_start_datetime):
                admindate = rad_pharm_seq[0][tag_Radiopharmaceutical_start_datetime].value
                adminDT = check_date_format(admindate)
            else:
                print("ERROR: PET pharm Admin date dicom tag is INCORRECT!")
            
            #scanDT = datetime.datetime.strptime(thedc[0x0009,0x100d].value,'%Y%m%d%H%M%S.%f')
            #adminDT = datetime.datetime.strptime(thedc[0x0009,0x103b].value,'%Y%m%d%H%M%S.%f')
            if ii == 0:
                adminDT0 = adminDT
                deltaT = scanDT - adminDT
            else:
                if adminDT != adminDT0:
                    deltaT = scanDT - adminDT0
                else:
                    deltaT = scanDT - adminDT

            # store dtSinceInj in unit of HOUR
            time_elapse_hr = deltaT.days*24. + (deltaT.seconds + deltaT.microseconds/1e6)/3600. # unit: seconds
            # print(scanDT, adminDT, adminDT0, deltaT, time_elapse_hr)
            self.PETdtSinceInj0.append(time_elapse_hr)
            tag_Units = Tag(0x0054, 0x1001)
            units_field = thedc[tag_Units]
            units = units_field.value
            print('image units:         {}'.format(units))

            # weight
            tag_Patient_Weight = Tag(0x0010, 0x1030)
            pat_weight_field = thedc[tag_Patient_Weight]
            pat_weight = pat_weight_field.value
            print('patient weight (kg): {}'.format(pat_weight))
            weight_units_factor = 1000.

            # injected dose
            tag_Radionuclide_Total_Dose = Tag(0x0018, 0x1074)
            tag_Radionuclide_Half_life = Tag(0x0018, 0x1075)
            try:
                rad_total_dose_field = rad_pharm_seq[0][tag_Radionuclide_Total_Dose]
            except:
                print('::O_O:: cannot get the total dose in the image DICOM header')
                return 0
            rad_total_dose = rad_total_dose_field.value
            rad_half_life = rad_pharm_seq[0][tag_Radionuclide_Half_life].value  # unit: seconds
            print('radionuclide total dose: {}'.format(rad_total_dose))

            # decay corrected the injected dose (or total dose)
            rad_lambda = mt.log(2) / rad_half_life  # unit: 1/sec
            decaycorr_rad_total_dose = rad_total_dose * mt.exp(-time_elapse_hr*3600. * rad_lambda)

            print('the decay-corrected total dose (MBq):    {}'.format(decaycorr_rad_total_dose))

            suv_multiplier = float(weight_units_factor * pat_weight) / float(decaycorr_rad_total_dose)
            self.PETSUV_factor0.append(suv_multiplier)
            print('suv_multiplier: {}'.format(suv_multiplier))

        # Set up Amide file names
        self.theAmideFnames = sorted(glob.glob('{}/VOIs_Amide/*.tsv'.format(self.PTdir)))
        self.theDay = sorted([int(re.findall(r'[\/_\w]+_Day([0-9]+)',ss)[0]) for ss in self.theAmideFnames])
        # self.PETdtSinceInj = sorted(self.PETdtSinceInj)
        # print(self.PETdtSinceInj0, self.PETSUV_factor0)
        self.PETdtSinceInj, self.PETSUV_factor = zip(*sorted(zip(self.PETdtSinceInj0, self.PETSUV_factor0)))


        #print self.theAmideFnames
        #print self.theDay
        # print(self.PETdtSinceInj)
        # print(self.PETSUV_factor)

    
    @staticmethod
    def RTFile2DF(fname):
        """
        - load a text file from Amide VOI output
        - get numerical output, variable titles, and ROI names
            
        """
        data = np.loadtxt(fname)
        with open(fname) as f:
            lines = f.readlines()

            tmpname = [re.findall('ROI:\t([\w\s()]+)\t',ll)[0] for ll in lines if re.findall('ROI:\t([\w()]+)',ll)]
            
            # make sure the organ names are in similar description as the src_organ name in the mysql database
            ROIname = [re.sub(r'T(\d)',r'Tumor\1',ss) for ss in tmpname]
            ROIname = map(lambda x: str.replace(x,'bladder','UrinaryBladder'),ROIname)
            
            for ll in lines:
                match = re.search(r'Frame',ll)
                # just find the first line with variable definition
                if match:
                    oh = re.split(r'\t+',ll)
                    # get rid of all white space or '#'
                    VARname = [re.sub(r'[#\s]+','',ss,flags=re.UNICODE) for ss in oh]
                    #print VARname
                    break
        
        df = pd.DataFrame(data,columns=VARname)

        # make sure all the organ are capitalized except for the all-capped string!
        CapROIname = [ss if ss.isupper() else upperfirst(ss) for ss in ROIname]
        df['OrganName'] = CapROIname
        df['Mean*Size(cm3)'] = df['Mean']*df['Size(mm^3)']/1000.
        return df

    def PMODFile2DF(self,fdir,ftag):
        # get the sub-directory of self.PMODdir
        allsubdir = os.walk(fdir).next()[1]

        # filter the directory name to 'TumorXXXX' only
        TumorName = [re.search(r'Tumor[0-9]+',dd).group() for dd in allsubdir if re.search(r'Tumor[0-9]+',dd)]

        self.theTumordf = pd.DataFrame()
        for tn in TumorName:
            allfiles = glob.glob('{}/{}/{}*.voistat'.format(fdir,tn,ftag))
            for ff in allfiles:
                tmpdf = pd.read_csv(ff,sep='\t',skiprows=1)
                nday = int(re.search(r'[\w\/]+{}_{}_Day([0-9]+)'.format(ftag,tn),ff).group(1))
                thesr = tmpdf.loc[1,['Total(AVR*VOL)','Averaged','Min','Max','Sd','Volume']]
                tmpstr = tmpdf.loc[0,['Total(AVR*VOL)','Averaged','Min','Max','Sd','Volume']]
                newindx = [s1 + s2 for s1,s2 in zip(thesr.index.tolist(),tmpstr.tolist())]
                thesr.index = newindx
                thesr['ImageDayNo.'] = nday
                thesr['PETdtSinceInjection_hr'] = self.PETdtSinceInj[self.theDay.index(nday)]
                thesr['PETSUV_factor'] = self.PETSUV_factor[self.theDay.index(nday)]
                thesr['RTFname'] = ff
                thesr['OrganName'] = tn
                self.theTumordf  = self.theTumordf.append(thesr,ignore_index=True)

        # convert numerica string to numeric float
        #self.theTumordf = self.theTumordf.convert_objects(convert_numeric=True)
        self.theTumordf = self.theTumordf.apply(pd.to_numeric,errors=('ignore'))
        #print self.theTumordf.dtypes
        self.theTumordf['Mean*Size(cm3)'] = self.theTumordf['Total(AVR*VOL)[(kBq/cc)*(ccm)]']*1000.


    def GetAllRTDF(self):
        RTFile2DF = ResTime.RTFile2DF
        
        combo = zip(self.theAmideFnames,self.PETdtSinceInj,self.theDay, self.PETSUV_factor)
        for ss,tt,dd, ff in combo:
            df_tmp = RTFile2DF(ss)
            df_tmp['RTFname'] = ss
            df_tmp['PETdtSinceInjection_hr'] = tt
            df_tmp['PETSUV_factor'] = ff
            df_tmp['ImageDayNo.'] = dd
            self.theRTdf = pd.concat([df_tmp,self.theRTdf],axis=0,ignore_index=True)

        # Combine RTdf with PMODDF
        if self.theTumordf is not None:
            self.theRTdf = self.theRTdf.append(self.theTumordf,ignore_index=True)
        #print df_tmp.count()
        #print self.theRTdf.count()


    def ComputeRT(self,isplot=False):
        organnameset = set(self.theRTdf['OrganName'].unique())
        self.theTBstr = organnameset.intersection(set(['Total','TotalBody','TB','WB'])).pop()
        theminDay = self.theRTdf['ImageDayNo.'].min()

        # decay correct total body PET counts
        df_totalbody = self.theRTdf[self.theRTdf['OrganName'] == self.theTBstr]
        # print df_totalbody.loc[df_totalbody['ImageDayNo.'] == theminDay,'Mean*Size(cm3)']
        dt = df_totalbody['PETdtSinceInjection_hr'].apply(lambda x: x - df_totalbody.loc[df_totalbody['ImageDayNo.'] == theminDay,'PETdtSinceInjection_hr'])
        dc_total = dt.iloc[:,0].apply(lambda x: df_totalbody.loc[df_totalbody['ImageDayNo.'] == theminDay,'Mean*Size(cm3)']*mt.exp(-x*self.img_isotope_lambda_p))
        # print(dc_total)
        
        self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr,'dc_dt'] = dt.iloc[:,0]
        # print self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr]
        self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr,'DecayCorr_Mean*Size(cm^3)'] = dc_total.iloc[:,0]
        # print self.theRTdf[self.theRTdf['OrganName'] == self.theTBstr]

        
        # combine VOI stats for SalivaryGlands
        if organnameset.intersection(set(['Salivary glands right','Salivary glands left'])):
            tmp1 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'Salivary glands right','Mean*Size(cm3)']
            tmp1.index = range(len(tmp1))
            tmp2 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'Salivary glands left','Mean*Size(cm3)']
            tmp2.index = range(len(tmp2))
            tmp3 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'Salivary glands right','PETdtSinceInjection_hr']
            tmp3.index = range(len(tmp3))
            tmp4 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'Salivary glands right','ImageDayNo.']
            tmp4.index = range(len(tmp4))
            tmp5 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'Salivary glands right', 'Mean']
            tmp5.index = range(len(tmp5))
            tmp6 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'Salivary glands left', 'Mean']
            tmp6.index = range(len(tmp6))
            list_series = [tmp1+tmp2,tmp5+tmp6, tmp3,tmp4]
            df_tmp = pd.concat(list_series,axis=1)
            df_tmp['OrganName'] = 'SalivaryGlands'
            # print(self.theRTdf.loc[self.theRTdf['OrganName'].isin(['Salivary glands right','Salivary glands left'])])
            self.theRTdf = self.theRTdf.append(df_tmp,ignore_index=True)
            # print(self.theRTdf.loc[self.theRTdf['OrganName'] == 'SalivaryGlands'])
        
        # combine VOI stats for Kidneys
        if organnameset.intersection(set(['RK','LK'])):
            tmp1 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'RK','Mean*Size(cm3)']
            tmp1.index = range(len(tmp1))
            tmp2 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'LK','Mean*Size(cm3)']
            tmp2.index = range(len(tmp2))
            tmp3 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'RK','PETdtSinceInjection_hr']
            tmp3.index = range(len(tmp3))
            tmp4 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'RK','ImageDayNo.']
            tmp4.index = range(len(tmp4))

            tmp5 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'RK', 'Mean']
            tmp5.index = range(len(tmp5))
            tmp6 = self.theRTdf.loc[self.theRTdf['OrganName'] == 'LK', 'Mean']
            tmp6.index = range(len(tmp6))

            list_series = [tmp1+tmp2,tmp5+tmp6, tmp3,tmp4]
            df_tmp = pd.concat(list_series,axis=1)
            df_tmp['OrganName'] = 'Kidney'
            self.theRTdf = self.theRTdf.append(df_tmp,ignore_index=True)

        # normalize all the organ PET signal by the decay-corrected total body count
        for dd in self.theDay:
            #print 'Day {}'.format(dd)
            oh = self.theRTdf[(self.theRTdf['OrganName'] == self.theTBstr ) & (self.theRTdf['ImageDayNo.'] == dd)]['DecayCorr_Mean*Size(cm^3)']
            self.theRTdf.loc[self.theRTdf['ImageDayNo.'] == dd,'P_normalizedPET'] = self.theRTdf.loc[self.theRTdf['ImageDayNo.'] == dd,'Mean*Size(cm3)'].apply(lambda x: 100.*x/oh.iloc[0])
            #print self.theRTdf[self.theRTdf['ImageDayNo.'] == dd]

        # make sure all normlaized PET signal is NOT negative
        self.theRTdf['P_normalizedPET'] = self.theRTdf['P_normalizedPET'].map(lambda x: x if x > 0.0 else 0.0)

        # compute SUV of the mean VOI signal
        self.theRTdf['PET_SUV'] = self.theRTdf.apply(PETSUV_func, axis=1)

        # exclude total body and salivary glands L & R
        self.theOrganName = [ss for ss in self.theRTdf['OrganName'].unique() if ss not in [self.theTBstr ,'Salivary glands right','Salivary glands left','RK','LK']]
        
        # set dtype == float32 even though there is a column of string. pandas will figure out what can be float32 and what cannot
        # need to make sure the data type in each column is set correctly or else groupby operation will not work properly since the column dtype is an 'object' instead of a float
        self.theResTimeHrDF = pd.DataFrame(columns = ['OrganName','a0','b0','a1','b1','ymax','r2','Residence Time (Bq-hr/Bq)','t0_hr','t1_hr','t2_hr','t3_hr','t4_hr','pInjAct0','pInjAct1','pInjAct2','pInjAct3','pInjAct4', 'SUV0','SUV1','SUV2','SUV3','SUV4', 'pIA_1_2TP_slope', 'pIA_2_3TP_slope', 'pIA_1_3TP_slope', 'SUV_1_2TP_slope','SUV_2_3TP_slope','SUV_1_3TP_slope'],index=range(len(self.theOrganName)+1),dtype='float32')

        for ii in range(len(self.theOrganName)):
            oname = self.theOrganName[ii]
            print 'organ name: {}'.format(oname)
            xdata = self.theRTdf.loc[self.theRTdf['OrganName'] == oname,'PETdtSinceInjection_hr'].as_matrix()
            ydata = self.theRTdf.loc[self.theRTdf['OrganName'] == oname,'P_normalizedPET'].as_matrix()
            sdata = self.theRTdf.loc[self.theRTdf['OrganName'] == oname, 'PET_SUV'].as_matrix()
            #print self.theRTdf.loc[self.theRTdf['OrganName'] == oname,'Mean*Size(cm3)'].as_matrix()

            # add data point (0.,0.)
            xdata = np.insert(xdata, 0, 0.0)
            ydata = np.insert(ydata, 0, 0.0)
            sdata = np.insert(sdata, 0, 0.0)

            #find the slope of the first 2 time points
            idx = np.argsort(xdata)
            pIA_2TP_slope = []
            SUV_2TP_slope = []
            for ii1,ii2 in self.TP_slope_idx:
                idx1 = np.where(idx == ii1)[0][0]
                idx2 = np.where(idx == ii2)[0][0]
                x1, x2 = xdata[idx1], xdata[idx2]
                y1, y2 = ydata[idx1], ydata[idx2]
                s1 ,s2 = sdata[idx1], sdata[idx2]
                pIA_2TP_slope.append((y2 - y1) / (x2 - x1))
                SUV_2TP_slope.append((s2 - s1) / (x2 - x1))

            # fit input: raw pixel dimension (not discrete 1,2,3,...) ==> pixdim = 1.0
            if isplot:
                fname = '{}/Summary/ResTimeFit_{}.pdf'.format(self.PTdir,self.theOrganName[ii])
                p1,r2,ymax,xfit,yfit = ivfit.FitBiExpo_ResTime(xdata,ydata,1.0,self.therapy_isotope_lambda_p,isplot=isplot,fname=fname)
            else:
                p1,r2,ymax,xfit,yfit = ivfit.FitBiExpo_ResTime(xdata,ydata,1.0,self.therapy_isotope_lambda_p)

            # data fit was done with normalized data (y/max); BE SURE TO MULTIPLE YMAX BACK TO THE ESTIMATED COEFFICENTS!
            rt = 0.01*ymax*(p1[0]/(p1[1] + self.therapy_isotope_lambda_p) + p1[2]/(p1[3] + self.therapy_isotope_lambda_p))
            dict_tmp = {'OrganName':oname,'a0':p1[0],'b0':p1[1],'a1':p1[2],'b1':p1[3],'ymax':ymax,'r2':r2,'Residence Time (Bq-hr/Bq)':rt}
            xdata_sort = xdata[idx]
            ydata_sort = ydata[idx]
            sdata_sort = sdata[idx]

            extra_keys = ['t{}_hr'.format(x) for x in range(len(xdata_sort))] + ['pInjAct{}'.format(x) for x in range(len(ydata_sort))] + \
                        ['SUV{}'.format(x) for x in range(len(sdata_sort))] + ['pIA_{}_{}TP_slope'.format(a,b) for a,b in self.TP_slope_idx] + \
                        ['SUV_{}_{}TP_slope'.format(a,b) for a,b in self.TP_slope_idx]
            extra_val = np.concatenate((xdata_sort, ydata_sort, sdata_sort, np.array(pIA_2TP_slope), np.array(SUV_2TP_slope)), axis=0).T
            dict_tac = dict(zip(extra_keys, extra_val))
            dict_tmp.update(dict_tac)
            self.theResTimeHrDF.iloc[ii,:] = pd.Series(dict_tmp)

        # calculate the residence time for whole-body (use a a*(1-exp(-b*x)  model fit)
        self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr,'PExcretion_Mean*Size(cm^3)'] = 100.*(self.theRTdf.loc[self.theRTdf['OrganName'] ==self.theTBstr,'DecayCorr_Mean*Size(cm^3)'] - self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr,'Mean*Size(cm3)'])/self.theRTdf.loc[self.theRTdf['OrganName'] ==self.theTBstr,'DecayCorr_Mean*Size(cm^3)']
        #print self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr,['PETdtSinceInjection_hr','Mean*Size(cm3)','DecayCorr_Mean*Size(mm^3)','PExcretion_Mean*Size(mm^3)']]
        xdata = self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr,'PETdtSinceInjection_hr'].as_matrix()
        ydata = self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr,'PExcretion_Mean*Size(cm^3)'].as_matrix()
        sdata = self.theRTdf.loc[self.theRTdf['OrganName'] == self.theTBstr, 'PET_SUV'].as_matrix()

        xdata = np.insert(xdata, 0, 0.0)
        ydata = np.insert(ydata, 0, 0.0)
        sdata = np.insert(sdata, 0, 0.0)

        # find the slope of the first 2 time points
        idx = np.argsort(xdata)
        pIA_2TP_slope = []
        SUV_2TP_slope = []
        for ii1, ii2 in self.TP_slope_idx:
            idx1 = np.where(idx == ii1)[0][0]
            idx2 = np.where(idx == ii2)[0][0]
            x1, x2 = xdata[idx1], xdata[idx2]
            y1, y2 = ydata[idx1], ydata[idx2]
            s1, s2 = sdata[idx1], sdata[idx2]
            pIA_2TP_slope.append((y2 - y1) / (x2 - x1))
            SUV_2TP_slope.append((s2 - s1) / (x2 - x1))

        if isplot:
            fname = '{}/Summary/ResTimeFit_TotalBody.pdf'.format(self.PTdir)
            p1,r2,ymax,xfit,yfit = ivfit.FitInvExpo(xdata,ydata,isplot=isplot,fname=fname)
        else:
            p1,r2,ymax,xfit,yfit = ivfit.FitInvExpo(xdata,ydata)
        rt = 0.01*((100.-ymax*p1[0])/self.therapy_isotope_lambda_p +  ymax*p1[0]/(p1[1] + self.therapy_isotope_lambda_p))
        if rt <  0.:
            print('::OH NO O_O:: Residence time of total body is NEGATIVE! Let\'s make some assumption and fit the model again!')
            # for the case where rt for total body was estimated to be neg (NO GOOD), usually due to missing time-point (fewer than 4 imaging time points)
            # assumption: add another data point at 120 hr, use the same %IA estimated for the last available time-point
            xdata = np.insert(xdata, len(xdata), 120.)
            ydata = np.insert(ydata, len(ydata), ydata[idx[-1]])
            sdata = np.insert(sdata, len(sdata), sdata[idx[-1]])
            if isplot:
                fname = '{}/Summary/ResTimeFit_TotalBody.pdf'.format(self.PTdir)
                p1, r2, ymax, xfit, yfit = ivfit.FitInvExpo(xdata, ydata, isplot=isplot, fname=fname)
            else:
                p1, r2, ymax, xfit, yfit = ivfit.FitInvExpo(xdata, ydata)
            rt = 0.01 * ((100. - ymax * p1[0]) / self.therapy_isotope_lambda_p + ymax * p1[0] / (p1[1] + self.therapy_isotope_lambda_p))
            idx = np.argsort(xdata) # re-sort the modified data

        dict_tmp = {'OrganName':'TotalBody','a0':p1[0],'b0':p1[1],'ymax':ymax,'r2':r2,'Residence Time (Bq-hr/Bq)':rt}
        xdata_sort = xdata[idx]
        ydata_sort = ydata[idx]
        sdata_sort = sdata[idx]

        extra_keys = ['t{}_hr'.format(x) for x in range(len(xdata_sort))] + ['pInjAct{}'.format(x) for x in range(len(ydata_sort))] + \
                     ['SUV{}'.format(x) for x in range(len(sdata_sort))] + ['pIA_{}_{}TP_slope'.format(a, b) for a, b in self.TP_slope_idx] + \
                     ['SUV_{}_{}TP_slope'.format(a, b) for a, b in self.TP_slope_idx]
        extra_val = np.concatenate((xdata_sort, ydata_sort, sdata_sort, np.array(pIA_2TP_slope), np.array(SUV_2TP_slope)), axis=0).T
        dict_tac = dict(zip(extra_keys, extra_val))
        dict_tmp.update(dict_tac)
        self.theResTimeHrDF.iloc[ii+1, :] = pd.Series(dict_tmp)

        # Keep the src organ name consistent: use the src name nomenclature defined in the srcname dct key in the .json file
        theTestlst = self.theResTimeHrDF.OrganName.tolist()
        theReflst = self.srcname_dict.keys()
        # themapdict = dict(zip(theTestlst, map(lambda p: Func_Match_Name1(p, theReflst), theTestlst)))
        # self.theResTimeHrDF['matched_OrganName'] = self.theResTimeHrDF['OrganName'].map(themapdict)
        self.theResTimeHrDF['matched_OrganName'] = self.theResTimeHrDF['OrganName'].apply(Func_Match_Name1, reflst=theReflst)
        self.theResTimeHrDF['matched_OrganName_lower'] = self.theResTimeHrDF['matched_OrganName'].str.lower()

        # # fill nan with other values for ResTimeHrDF
        # self.theResTimeHrDF.fillna(0.0, inplace=True)

    def ComputeOrganDose(self,theGeoIdlst,theSimPkg,InjDoseMBq):
        # Pull a list of S-value from the database
        qr = "SELECT geo_id, b.src_organ, a.target_organ, a.SV_mean from DoseInfo a JOIN SimInfo b ON a.sim_id = b.sim_id WHERE geo_id in ({}) and simpkg = '{}'".format(','.join('"{0}"'.format(w) for w in theGeoIdlst),theSimPkg)
        df_all = pd.read_sql(qr,self.mysqlcon)

        # # Keep the src organ name consistent: use the src name nomenclature defined in the srcname dct key in the .json file
        # theTestlst = self.theResTimeHrDF.OrganName.tolist()
        # theReflst = df_all.src_organ.unique().tolist()
        # themapdict = dict(zip(theTestlst,map(lambda p: NameMatchFunc(p,theReflst),theTestlst)))
        # self.theResTimeHrDF['matched_OrganName'] = self.theResTimeHrDF['OrganName'].map(themapdict)
        # self.theResTimeHrDF['matched_OrganName_lower'] = self.theResTimeHrDF['matched_OrganName'].str.lower()

        # make sure the organ name are all lower case
        df_all['src_organ_lower'] = df_all['src_organ'].str.lower()
        df_join1 = pd.merge(self.theResTimeHrDF,df_all,left_on = 'matched_OrganName_lower',right_on = 'src_organ_lower',how='left')
        #print len(df_join1.src_organ.unique()), len(theTestlst)


        # create another df without 'totalbody' as src_organ
        df_dose = df_join1[df_join1['src_organ'] != 'TotalBody']

        # take heed for computing organ dose from 'Remainder body' as a source organ
        # get organ mass for a given geometry
        qr = "SELECT geo_id, OrganName,Mass_g,Volume_cm3 from GeoInfo WHERE geo_id in ({})".format(','.join('"{0}"'.format(w) for w in theGeoIdlst))
        df_mass = pd.read_sql(qr,self.mysqlcon)
        
        # find mass and volume of the remaining body
        # date frame series iget(i) returns the i-th value in the series by location
        solst = df_join1.loc[df_join1['src_organ'] != 'TotalBody','src_organ'].unique()
        mass_RB = (df_mass.loc[df_mass['OrganName'] == 'TotalBody','Mass_g'] - df_mass[df_mass['OrganName'].isin(solst)]['Mass_g'].sum()).iloc[0]
        vol_RB = (df_mass.loc[df_mass['OrganName'] == 'TotalBody','Volume_cm3'] - df_mass[df_mass['OrganName'].isin(solst)]['Volume_cm3'].sum()).iloc[0]
        #print mass_RB, vol_RB

        # join 2 tables
        df_join2 = pd.merge(df_all,df_mass,left_on = 'src_organ',right_on = 'OrganName',how='left')
        df_join2['SVxMass'] = df_join2['SV_mean']*df_join2['Mass_g']

        # Estimate the SV for 'Remainder body': "Fundamentals of Nuclear Medicine Dosimetry" by M. Stabin, p104
        # R. Cloutier et al., JNM 1973 "Calculating the radiation dose to an organ", see eq 14
        # SV_RB = SV(rk<-rTB)*M_TB/M_RB - sum_over_h[SV(rk<-rh)*M_h/M_RB]
        sr_TB = df_join2.query('src_organ == "TotalBody"').groupby('target_organ').agg({'SVxMass':np.mean})
        sr_others = df_join2.query('src_organ != "TotalBody"').groupby('target_organ').agg({'SVxMass':np.sum})
        df_RB = sr_TB.subtract(sr_others)/mass_RB
        df_RB.columns = ['SV_mean']

        # new df with labeling info in the columns under index names and append the RB dose dataframe with dose info of the other src organs
        df_RB.reset_index(level=0,inplace=True)
        ResTime_RB = self.theResTimeHrDF.loc[self.theResTimeHrDF['OrganName']=='TotalBody','Residence Time (Bq-hr/Bq)'].iloc[0] - self.theResTimeHrDF.loc[self.theResTimeHrDF['OrganName']!='TotalBody','Residence Time (Bq-hr/Bq)'].sum()
        df_RB['Residence Time (Bq-hr/Bq)'] = ResTime_RB
        df_RB['src_organ'] = 'RemainderBody'
        df_dose = df_dose.append(df_RB,ignore_index=True)

        # compuate the organ dose: multiply the S-value to residence time and initial injection dose
        # S-value: Svalue (mGy/MBq-s) * ResTime(MBq-hr/MBq) * A0(MBq)
        df_dose['OrganDose(mGy)'] = df_dose['SV_mean']*df_dose['Residence Time (Bq-hr/Bq)']*3600*InjDoseMBq  # unit: mGy/MBq * (MBq) = mGy

        # calculate organ dose from all SOURCE ORGANS except 'TotalBody'
        self.theOrganDosemGy = df_dose.groupby('target_organ').agg({'OrganDose(mGy)':np.sum})
        self.theOrganDosemGy.reset_index(level=0,inplace=True)
        #print self.theOrganDosemGy

        # append via data frame not series (http://stackoverflow.com/questions/24284342/insert-a-row-to-pandas-dataframe)
        self.theResTimeMass = df_join1.groupby('src_organ').agg({'Residence Time (Bq-hr/Bq)':np.mean})
        self.theResTimeMass['Mass_g'] = df_join2.groupby('src_organ').agg({'Mass_g':np.mean})
        self.theResTimeMass['Volume_cm3'] = df_join2.groupby('src_organ').agg({'Volume_cm3':np.mean})
        self.theResTimeMass = self.theResTimeMass.append(pd.DataFrame({'Residence Time (Bq-hr/Bq)':ResTime_RB,'Mass_g':mass_RB,'Volume_cm3':vol_RB},index=['RemainderBody']))
        self.theResTimeMass.reset_index(level=0,inplace=True)

    def ComputeEffectiveDose(self, wt_type=''):
        fname = '{}/OrganData/EffectiveDoseWt.csv'.format(self.Frdatadir)
        self.df_EDWt = pd.read_csv(fname, header=0)
        EDwt_reflst = self.df_EDWt['Target organ'].tolist()

        # match organ name with the wt factor organ name
        if self.theOrganDosemGy is not None: # make sure organ dose DF is defined and populated
            self.theOrganDosemGy['Target_name_match'] = self.theOrganDosemGy['target_organ'].apply(Func_Match_Name2, reflst=EDwt_reflst)
            tmp_df = self.theOrganDosemGy.merge(self.df_EDWt, how='inner',left_on='Target_name_match', right_on='Target organ')
            tmp_df['ED*dose'] = tmp_df.apply(lambda row: row['OrganDose(mGy)']*row['ICRP 103 wt'],axis=1)
            self.theED = tmp_df['ED*dose'].sum()/1000. #unit: Sv
        else:
            print('No organ dose was calculated!! O_O')



