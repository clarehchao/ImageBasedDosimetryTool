from __future__ import division
import re
import pandas as ps
import numpy as np
import os

"""
Developer Note:
09.23.2015: as you notice, create table, insert into tabls are not table specific. will be nice to write a separate class defining how each table is inserted and setup etc
            one can pass in each 'table' specific class to work with DoseDB toolset in a more general way without having to define the same function for different tables
"""

def get_DB_auth_info(db_auth_dir, db_user_name=''):
    db_pw_path = os.path.join(db_auth_dir, db_user_name + '.auth')
    f = open(db_pw_path)
    db_pw = f.readline().strip()
    f.close()
    return db_pw


def ReadFile2DF(fname,isheader):
    if isheader:
        df = ps.read_table(fname,index_col=0)  # index_col = specify the column to use as the row labels of the data frame
    else:
        df = ps.read_table(fname,header=None,index_col=0)  # header = row number to use as the column names, header = None if the file has no header row  
    return df
    
def CreateTableDB_DoseSimInfo(con):
    with con: # 'with' keyword automatically release the resource (i.e. close the db, catch error)
        cur = con.cursor()
        # create table for simulation info of the dose table
        cur.execute("DROP TABLE IF EXISTS SimInfo")
        cur.execute("CREATE TABLE SimInfo(sim_id INT UNSIGNED NOT NULL PRIMARY KEY,simpkg VARCHAR(25) NOT NULL,geo_id VARCHAR(100) NOT NULL, \
        src_particle VARCHAR(25) NOT NULL, src_organ VARCHAR(25) NOT NULL,Nevent INT NULL, Nrun INT NULL);")
                     

        # create table for dose
        cur.execute("DROP TABLE IF EXISTS DoseInfo")
        cur.execute("CREATE TABLE DoseInfo(sim_id INT UNSIGNED NOT NULL,target_organ VARCHAR(25) NOT NULL,SV_mean DECIMAL(20,15) NULL, SV_std DECIMAL(20,15) NULL, \
        PRIMARY KEY(sim_id,target_organ));")
                     
                     
def CreateTableDB_GeoInfo(con):
    with con: # 'with' keyword automatically release the resource (i.e. close the db, catch error)
        cur = con.cursor()
        # create table for geometry information e.g. mass, volume
        cur.execute("DROP TABLE IF EXISTS GeoInfo")
        cur.execute("CREATE TABLE GeoInfo(id INT NOT NULL AUTO_INCREMENT,geo_id VARCHAR(100) NOT NULL, OrganName VARCHAR(100) NOT NULL, \
        Volume_cm3 DOUBLE NULL,Mass_g DOUBLE NULL,PRIMARY KEY(id));")

def CreateTableDB_MIBGPTInfo(con):
    with con: # 'with' keyword automatically release the resource (i.e. close the db, catch error)
        cur = con.cursor()
        # create table for geometry information e.g. mass, volume
        cur.execute("DROP TABLE IF EXISTS MIBGPTInfo")
        cur.execute("CREATE TABLE MIBGPTInfo(id INT NOT NULL AUTO_INCREMENT, pt_id VARCHAR(100) NOT NULL,PETCTDir VARCHAR(500) NOT NULL, ResTimeDir VARCHAR(500), \
        geo_id VARCHAR(500) NOT NULL,I131MIBGDose_mCi DOUBLE NULL,PRIMARY KEY(id));")

def CreateTableDB_ResTimeInfo(con):
    with con: # 'with' keyword automatically release the resource (i.e. close the db, catch error)
        cur = con.cursor()
        # create table for geometry information e.g. mass, volume
        cur.execute("DROP TABLE IF EXISTS ResTimeInfo")
        cur.execute("CREATE TABLE ResTimeInfo(id INT NOT NULL AUTO_INCREMENT,pt_id VARCHAR(100) NOT NULL,OrganName VARCHAR(500) NOT NULL, \
        a0 DOUBLE NULL,b0 DOUBLE NULL,a1 DOUBLE NULL,b1 DOUBLE NULL,ymax DOUBLE NULL,r2 DOUBLE NULL, ResTime_BqhrPerBq DOUBLE NULL, \
        t0_hr DOUBLE NULL, t1_hr DOUBLE NULL,t2_hr DOUBLE NULL,t3_hr DOUBLE NULL,t4_hr DOUBLE NULL, pInjAct0 DOUBLE NULL, pInjAct1 DOUBLE NULL, \
        pInjAct2 DOUBLE NULL, pInjAct3 DOUBLE NULL, pInjAct4 DOUBLE NULL, SUV0 DOUBLE NULL, SUV1 DOUBLE NULL, SUV2 DOUBLE NULL, SUV3 DOUBLE NULL, \
        SUV4 DOUBLE NULL, pIA_1_2TP_slope DOUBLE NULL, pIA_2_3TP_slope DOUBLE NULL, pIA_1_3TP_slope DOUBLE NULL, SUV_1_2TP_slope DOUBLE NULL, \
        SUV_2_3TP_slope DOUBLE NULL, SUV_1_3TP_slope DOUBLE NULL, PRIMARY KEY(id));")

def CreateTableDB_AbsorbedDoseInfo(con):
    with con: # 'with' keyword automatically release the resource (i.e. close the db, catch error)
        cur = con.cursor()
        # create table for geometry information e.g. mass, volume
        cur.execute("DROP TABLE IF EXISTS AbsorbedDoseInfo")
        cur.execute("CREATE TABLE AbsorbedDoseInfo(id INT NOT NULL AUTO_INCREMENT, pt_id VARCHAR(100) NOT NULL,TargetOrgan VARCHAR(500) NOT NULL,AbsorbedDose_mGy DOUBLE NULL, \
        PRIMARY KEY(id));")

def CreateTableDB_EDInfo(con):
    with con: # 'with' keyword automatically release the resource (i.e. close the db, catch error)
        cur = con.cursor()
        # create table for geometry information e.g. mass, volume
        cur.execute("DROP TABLE IF EXISTS EDInfo")
        cur.execute("CREATE TABLE EDInfo(id INT NOT NULL AUTO_INCREMENT, pt_id VARCHAR(100) NOT NULL, geo_id VARCHAR(100) NOT NULL, \
        EffectiveDose_Sv DOUBLE NULL, PRIMARY KEY(id));")
                     
def Insert2DB_GeoInfo(con,df,vargs):
    geo_id = vargs
    
    df_fill = df.fillna(0)
    with con:
        cur = con.cursor()
        tmpdf = df_fill
        tmpdf['geo_id'] = geo_id
        tmplst = list(tmpdf.loc[:,['geo_id','OrganName','Volume (cm3)','Mass (g)']].itertuples(index=False))
        cur.executemany("INSERT INTO GeoInfo(geo_id,OrganName,Volume_cm3,Mass_g) VALUES(%s,%s,%s,%s);",tmplst)
        
def Insert2DB_MIBGPTInfo(con,df):
    # make sure there is no NAN in the dataframe
    df_fill = df.fillna(0)
    with con:
        cur = con.cursor()
        tmplst = list(df_fill.loc[:,['PTid','PETCTdir','ResTimedir','geo_id','I131MIBGDose_mCi']].itertuples(index=False))
        cur.executemany("INSERT INTO MIBGPTInfo(pt_id,PETCTDir,ResTimeDir,geo_id,I131MIBGDose_mCi) VALUES(%s,%s,%s,%s,%s);",tmplst)

def Insert2DB_ResTimeInfo(con,df,vargs):
    ptid = vargs

    # make sure there is no NAN in the dataframe
    df_fill = df.fillna(0.0)
    with con:
        cur = con.cursor()
        tmpdf = df_fill
        tmpdf['ptid'] = ptid
        TP_slope_idx_lst = [(1,2),(2,3),(1,3)]
        colname_lst = ['ptid','matched_OrganName','a0','b0','a1','b1','ymax','r2','Residence Time (Bq-hr/Bq)'] + ['t{}_hr'.format(x) for x in range(5)] + \
                      ['pInjAct{}'.format(x) for x in range(5)] + ['SUV{}'.format(x) for x in range(5)] + \
                      ['pIA_{}_{}TP_slope'.format(a,b) for a,b in TP_slope_idx_lst] + ['SUV_{}_{}TP_slope'.format(a,b) for a,b in TP_slope_idx_lst]
        tmplst = list(tmpdf.loc[:, colname_lst].itertuples(index=False))
        cur.executemany("INSERT INTO ResTimeInfo(pt_id,OrganName,a0,b0,a1,b1,ymax,r2,ResTime_BqhrPerBq,t0_hr,t1_hr,t2_hr,t3_hr,t4_hr,pInjAct0,pInjAct1,pInjAct2,pInjAct3,pInjAct4,SUV0,SUV1,SUV2,SUV3,SUV4,pIA_1_2TP_slope,pIA_2_3TP_slope,pIA_1_3TP_slope, SUV_1_2TP_slope, SUV_2_3TP_slope,SUV_1_3TP_slope) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);",tmplst)

def Insert2DB_AbsorbedDoseInfo(con,df,vargs):
    ptid = vargs

    # make sure there is no NAN in the dataframe
    df_fill = df.fillna(0)
    with con:
        cur = con.cursor()
        tmpdf = df_fill
        tmpdf['ptid'] = ptid
        tmplst = list(tmpdf.loc[:,['ptid','target_organ','OrganDose(mGy)']].itertuples(index=False))
        cur.executemany("INSERT INTO AbsorbedDoseInfo(pt_id,TargetOrgan,AbsorbedDose_mGy) VALUES(%s,%s,%s);",tmplst)

    
def Insert2DB_DoseSimInfo(con,df,varags):
    target_organs = df.index
    src_organs = df.columns
    simpkg,geo_id,src_particle,last_sim_id = varags
    with con:
        cur = con.cursor()
        for ii in range(len(src_organs)):
            tmp1 = [last_sim_id+ii+1,simpkg,geo_id,src_particle,src_organs[ii]]
            cur.execute("INSERT INTO SimInfo(sim_id,simpkg,geo_id,src_particle,src_organ) VALUES(%s,%s,%s,%s,%s);",tmp1)
            for jj in range(len(target_organs)):
                tmp2 = [last_sim_id+ii+1,target_organs[jj],df[src_organs[ii]][jj]]
                cur.execute("INSERT INTO DoseInfo(sim_id,target_organ,SV_mean) VALUES(%s,%s,%s);",tmp2)
                
def Insert2DB_fancydf_DoseSimInfo(con,df,varags):
    target_organs = df.index
    last_sim_id = ReadDBSize(con)
    start_sim_id = int(last_sim_id) + 1
    #simpkg,geo_id,src_particle,src_organ,Nevent,Nrun = varags
    with con:
        cur = con.cursor()
        varags.insert(0,start_sim_id)
        cur.execute("INSERT INTO SimInfo(sim_id,simpkg,geo_id,src_particle,src_organ,Nevent,Nrun) VALUES(%s,%s,%s,%s,%s,%s,%s);",varags)
        for i in range(len(target_organs)):
            to = target_organs[i]
            tmp = [start_sim_id,to,df[1][to],df[2][to]]
            cur.execute("INSERT INTO DoseInfo(sim_id,target_organ,SV_mean,SV_std) VALUES(%s,%s,%s,%s);",tmp)
            
def Insert2DB_fancydict_DoseSimInfo(con,dct,varags):
    last_sim_id = ReadDBSize(con)
    start_sim_id = int(last_sim_id) + 1
    #simpkg,geo_id,src_particle,src_organ,Nevent,Nrun = varags
    with con:
        cur = con.cursor()
        varags.insert(0,start_sim_id)
        #print varags
        cur.execute("INSERT INTO SimInfo(sim_id,simpkg,geo_id,src_particle,src_organ,Nevent,Nrun,G4AppName,G4SimDataDir) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s);",varags)
        for to,val in dct.items(): 
            tmp = [start_sim_id,to] + val
            #print tmp
            cur.execute("INSERT INTO DoseInfo(sim_id,target_organ,SV_mean,SV_std) VALUES(%s,%s,%s,%s);",tmp)


def Insert2DB_EDInfo(con, vargs):
    #pt_id, geo_tag, ED = vargs
    with con:
        cur = con.cursor()
        cur.execute("INSERT INTO EDInfo(pt_id,geo_id,EffectiveDose_Sv) VALUES(%s,%s,%s);", vargs)
       
def ReadDBSize(con):
    with con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM SimInfo")
        tmp = cur.fetchone()
    return tmp[0]
    
def GetNevent(fname):
    # Determine the number of total events in the dose log file
    with open(fname,'r') as f:
        nevent = re.findall(r'Number of Events in this run: ([\w]+)',f.read())[0]
    return int(nevent)

def Log2Data(dosedir,runlist):
    neventlist = []
    for r in runlist:
        logfile = '{}/Run{}/log.txt'.format(dosedir,r)
        neventlist.append(GetNevent(logfile))
    return neventlist
    
def FindIntermIrun(dosedir,run1,inc1,neventlist):
    for i in range(1,11):
        irun = run1 + i*inc1
        logfile = '{}/Run{}/log.txt'.format(dosedir,irun)
        neventeh = GetNevent(logfile)
        if neventeh == neventlist[1]:
            break
    return irun
     
def Fname2Data(fname,simfdir):
    # example: SVstatsufh10f_1_ufh10f_1_UBCont_I131_Run1-69_rinc20_4_G4.9.6.p02.txt
    tmp1 = re.findall(r'SVstats([\w-]+)_(\w+)_(\w+)_Run(\d+)-(\d+)_rinc([\w-]*)_([\w.]+).txt',fname)  #this separate the filename into sim info
    print tmp1
    junk,src_organ,src_particle,run1,run2,inc,pkg = tmp1[0]
    tmp2 = re.findall(r'(.+?)_\1+',fname)  #this finds the repeated str, e.g. \1 means 'repeat the pattern in the space before \1, which is '(.+?)_'
    geo_id = tmp2[0]
    
    # get Nevent info
    ddir = '{}/GEO_{}/SRCMP_{}/{}/{}'.format(simfdir,geo_id,geo_id,src_organ,src_particle)
    run_l = [int(run1),int(run2)]
    nevent_l = Log2Data(ddir,run_l)
    #print nevent_l
    
    # check for multiple entris of increment
    check = re.findall(r'(\d+)_(\d+)',inc)
    if check:  # found something in check
        inc_l = [int(s) for s in check[0]]
        iirun = FindIntermIrun(ddir,run_l[0],inc_l[0],nevent_l)
        run1_ar = np.array([run_l[0],iirun])
        run2_ar = np.array([iirun-run_l[0],run_l[1]+inc_l[1]-1])
        y = run2_ar - run1_ar + 1
        #print iirun,run1_ar,run2_ar,y
        Nrun = int(sum(y/np.array(inc_l)))  # number of simulations
        Nevent = int(inc_l[0]*nevent_l[0])  # number of particles simulated PER event
        #np.array(inc_l)*np.array(nevent_l)
    else:
        iinc = int(inc)
        # original: Nrun = (run_l[1] + iinc - 1 - run_l[0] + 1)/iinc, simplified to the below
        Nrun = int((run_l[1] + iinc - run_l[0])/iinc)
        Nevent = int(nevent_l[0]*iinc)
    
    return [pkg,geo_id,src_particle,src_organ,Nevent,Nrun]
    
def CheckTableExist(con,tablename):
    cur = con.cursor()
    mysqlstm = "SHOW TABLES LIKE '{}'".format(tablename)
    result = cur.execute(mysqlstm)
    if result:
        return True
    else:
        return False


