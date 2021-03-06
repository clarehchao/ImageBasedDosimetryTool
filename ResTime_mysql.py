#! /usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as ps
import ResTime as rt
from ImVolTool import FileIO as fi
from ImVolTool import xyz2vol as xv
import sys
from DBTool import DoseDB as ddb
import MySQLdb as mdb
import matplotlib.pyplot as plt
import seaborn as sns



if __name__ == '__main__':
    finput = sys.argv[1]
    param_dict = fi.LoadInputParameter(finput)


    # set up the MIBG patient file directory
    ptdir = '{}/PT{:0>4d}'.format(param_dict['NBptdir'],param_dict['PTid'])
    thePETdir = '{}/IM'.format(ptdir)
    theAMIDEVOIdir = '{}/VOIs_Amide'.format(ptdir)
    I131MIBGinjDoseMBq = param_dict['I131MIBGinjDosemCi']*37.

    # Initialize the MySQL connection
    # this works on Higgs as if one specifies 'localhost', it tries to connect via default /tmp/mysql.socket
    # see http://stackoverflow.com/questions/4662364/cant-connect-to-localhost-using-pythons-mysqldb
    #con = mdb.connect('localhost','testuser','test000','UCSFDoseDB')
    # con = mdb.connect('127.0.0.1','testuser','test000','UCSFDoseDB')
    db_pw = ddb.get_DB_auth_info(param_dict['DB_auth_dir'],param_dict['DB_usr'])
    con = mdb.connect(host='127.0.0.1', user=param_dict['DB_usr'],passwd=db_pw, db='UCSFDoseDB')

    # create a ResTime Constructor with initial variable setup
    theRT = rt.ResTime(ptdir,param_dict['petdirs'],param_dict['therapy_isotope_halflife_day']*24.,
                       param_dict['img_isotope_halflife_day']*24., param_dict['srcname'], param_dict['VHDMSDdir'],
                       param_dict['DB_auth_dir'], param_dict['DB_usr'])
    theRT.InitialSetUp()
    #print theRT.PETdtSinceInj

    # Read PMOD files if there are PMOD files contour available for the patient
    if param_dict['pmodftag']:
        print 'PMOD files exist! set up the pmod file directory!'
        pmoddir = '{}/PMODData'.format(ptdir)
        theRT.PMODFile2DF(pmoddir,param_dict['pmodftag'])

    # make xxxx/Summary directoy for saving plots
    if param_dict['isFitPlot']:
        fdir = '{}/Summary'.format(ptdir)
        xv.MakeDir(fdir)

    # Combine all the PET time point data for all source organs
    theRT.GetAllRTDF()
    #print theRT.theRTdf

    # Compute the residence for each source organ (get time activity curve and Bi-expo fit)
    theRT.ComputeRT(isplot=param_dict['isFitPlot'])
    # print(theRT.theResTimeHrDF)
    
    
    if not ddb.CheckTableExist(con,'ResTimeInfo'):
        ddb.CreateTableDB_ResTimeInfo(con)
    ddb.Insert2DB_ResTimeInfo(con,theRT.theResTimeHrDF,param_dict['PTid'])

    # # fname = '{}/theRTdf_PT{}.csv'.format(theAMIDEVOIdir, param_dict['PTid'])
    # # theRT.theRTdf.to_csv(fname)

    # compute the final organ dose in mGy considering the MIBG therapy dose and dose contribution from all source organs
    I131MIBGinjDoseMBq = param_dict['I131MIBGinjDosemCi']*37.
    theRT.ComputeOrganDose([param_dict['geotag']],param_dict['G4SimPkg'],I131MIBGinjDoseMBq)
    # print(theRT.theOrganDosemGy)
    # fname = '{}/Summary/theOrganDosemGy_PT{}.csv'.format(ptdir, param_dict['PTid'])
    # theRT.theOrganDosemGy.to_csv(fname)
    
    if not ddb.CheckTableExist(con,'AbsorbedDoseInfo'):
        ddb.CreateTableDB_AbsorbedDoseInfo(con)
    ddb.Insert2DB_AbsorbedDoseInfo(con,theRT.theOrganDosemGy,param_dict['PTid'])

    theRT.ComputeEffectiveDose()
    # print(theRT.theED)

    if not ddb.CheckTableExist(con,'EDInfo'):
        ddb.CreateTableDB_EDInfo(con)
    ddb.Insert2DB_EDInfo(con,[param_dict['PTid'], param_dict['geotag'], theRT.theED])

    if not ddb.CheckTableExist(con,'MIBGPTInfo'):
        print 'table does not exist!'
        ddb.CreateTableDB_MIBGPTInfo(con)
    ptinfo_df = ps.DataFrame([{'PTid':param_dict['PTid'],'PETCTdir':thePETdir,'ResTimedir':theAMIDEVOIdir,'geo_id':param_dict['geotag'],'I131MIBGDose_mCi':param_dict['I131MIBGinjDosemCi']}])
    #ptinfo_df = ptinfo_df.convert_objects(convert_numeric=True)
    ddb.Insert2DB_MIBGPTInfo(con,ptinfo_df)
    
    if param_dict['isDosePlot']:
         # plot the target organ vs absorbed dose
        df_OD = theRT.theOrganDosemGy
        df_OD['OrganDoseGy'] = df_OD['OrganDose(mGy)'].apply(lambda x: x*0.001)
        sns.set_context('talk')
        g = sns.factorplot(x='target_organ',y='OrganDoseGy',data=df_OD.sort_values(by='OrganDoseGy'),kind='bar',color=sns.color_palette('GnBu_d')[5],size=9,aspect=1.5,edgecolor='1.0')
        ax = g.fig.get_axes()[0]
        #ax.set_xticklabels(df_OD['TargetOrgan'].tolist(),rotation=20,fontsize=10)
        ax.set_xlabel('')
        ax.set_xticklabels(df_OD.sort_values(by='OrganDoseGy').target_organ.tolist(),rotation=20,fontsize=10)
        ax.set_ylabel('Absorbed Dose (Gy)',fontsize=20)
        g.despine(left=True)
        #ax.tick_params(axis='x',which='major',pad=-50)
        fname = '{}/Summary/BarPlot_ADvsTargetOrgan.pdf'.format(ptdir)
        plt.savefig(fname)
        #plt.show()
    #     # plot residence time and organ mass
        df_RT = theRT.theResTimeMass
        thedata = df_RT[~df_RT['src_organ'].isin(['RemainderBody','TotalBody'])]
        #print df_RT[df_RT['src_organ'] == 'RemainderBody']
        sns.set(context="talk")
        f, (ax1,ax2) = plt.subplots(2,1,figsize=(11,9),sharex=True)
        sns.barplot(x='src_organ',y='Mass_g',data=thedata,color=sns.color_palette('Set3')[3],ax=ax1,edgecolor='1.0')
        sns.barplot(x='src_organ',y='Residence Time (Bq-hr/Bq)',data=thedata,color=sns.color_palette('Set3')[4],ax=ax2,edgecolor='1.0')
        ax1.set_ylabel('Mass (g)',fontsize=15)
        ax1.set_xlabel('')
        ax1.set_xticklabels(thedata['src_organ'].tolist(),rotation=30,fontsize=10)
        ax2.set_ylabel('Residence Time (Bq-hr/Bq)',fontsize=15)
        ax2.set_xlabel('')
        ax2.set_xticklabels(thedata['src_organ'].tolist(),rotation=30,fontsize=10)
        sns.despine(bottom=True)
        #plt.setp(f.axes, yticks=[])
        plt.tight_layout(h_pad=3)
        fname = '{}/Summary/BarPlot_MassANDRTvsOrgan.pdf'.format(ptdir)
        plt.savefig(fname)
        #plt.show()
        plt.close('all')
