#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 3/26/18

@author: shuang
@goal: analyze the relationship between the two-point slope and the actual residence time

"""

import pandas as pd
import MySQLdb as mdb
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np



def query_mysql(con,qr):
    """
    :param con: connection
    :param qr: mysql query
    :return: a pandas dataframe of the data from the query
    """
    df = pd.read_sql(qr,con)
    return df

def srcname_redef(s):
    if s.lower().find('adrenal') != -1:
        return 'Adrenal Glands'
    elif s.lower().find('salivary') != -1:
        return 'Salivary Glands'
    elif s.lower().find('tumor') != -1:
        return 'Tumor'
    elif s.lower().find('bladder') != -1:
        return 'Urinary Bladder'
    else:
        return s

data_dir = '/Users/shuang/Documents/Proj_Neuroblastoma/Data/ResTime_Analysis'

# con = mdb.connect('localhost','testuser','test000','UCSFDoseDB')
con = mdb.connect(host='127.0.0.1', user='root', passwd='TWvachian81', db='UCSFDoseDB')

# get the residence time info
qr_rt1 = "SELECT * FROM ResTimeInfo"
df1 = query_mysql(con, qr_rt1)
df1['OrganName_mod'] = df1['OrganName'].map(srcname_redef)

# get the absorbed dose info
qr_rt2 = "SELECT * FROM AbsorbedDoseInfo"
df2 = query_mysql(con, qr_rt2)
df2['OrganName_mod'] = df2['TargetOrgan'].map(srcname_redef)
df2['AbsorbedDose_Gy'] = df2['AbsorbedDose_mGy'].map(lambda x: x/1000.)


# # save the raw data to csv
# csv_fname = '{}/Organ_Tumor_ResTime_2TPslope_all.csv'.format(data_dir)
# df1.drop(['id','OrganName_mod'], axis=1).to_csv(csv_fname, index=False)


# selet the df where it's not tumor and plot ResTime vs Organ
df3 = df1[(df1['OrganName_mod'] != 'Tumor') & (df1['OrganName_mod'] != 'TotalBody')]

# TODO: need to fix the plot and make sure it's showing things clearly!
# TODO: PT10 has a total body RT with neg value, .... how do i remedy this? the patient only had 3 imaging time point

sns.set(font='Arial')
sns.set_context('paper')

plt.figure(figsize=(11,8.5))
# f, ax = plt.subplots()
# f.set_size_inches(11,8.5)
g1 = sns.stripplot(x='ResTime_BqhrPerBq', y='OrganName_mod', data=df3, dodge=True, jitter=True, alpha=0.5, zorder=1, palette='dark')
g2 = sns.pointplot(x='ResTime_BqhrPerBq', y='OrganName_mod', data=df3, dodge=0.532, join=False, palette='dark', markers='d', scale=1.2, ci=None, estimator=np.mean)
# g1.set_xlim(0,11.)
g1.set_xlabel('Residence time (Bq*hr/Bq)',{'weight': 'bold', 'size': 15})
g1.set_ylabel('Organ', {'weight': 'bold', 'size': 15})
# print(df3.ix[df3['OrganName_mod'] == 'Liver','ResTime_BqhrPerBq'].mean())

# stackoverflow post on get_xticklabel being emtpy string: https://stackoverflow.com/questions/11244514/modify-tick-label-text
xtickll = g1.get_xticks().tolist()

g1.set_xticklabels(xtickll, {'weight': 'bold', 'size': 13})
g1.set_yticklabels(g1.get_yticklabels(), {'weight': 'bold', 'size': 13})
plt.tight_layout()
fig_name = '{}/ResTime_Organs.pdf'.format(data_dir)
plt.savefig(fig_name)
# plt.show()


# selet the df where it's not tumor and plot Absorbed dose vs Organ
df4 = df2[(df2['OrganName_mod'] != 'Tumor') & (df2['OrganName_mod'] != 'TotalBody')]

sns.set(font='Arial')
sns.set_context('paper')

plt.figure(figsize=(11,8.5))
# f, ax = plt.subplots()
# f.set_size_inches(11,8.5)
g1 = sns.stripplot(x='AbsorbedDose_Gy', y='OrganName_mod', data=df4, dodge=True, jitter=True, alpha=0.5, zorder=1, palette='dark')
g2 = sns.pointplot(x='AbsorbedDose_Gy', y='OrganName_mod', data=df4, dodge=0.532, join=False, palette='dark', markers='d', scale=1.2, ci=None)
# g1.set_xlim(0,5)
g1.set_xlabel('Absorbed Dose (Gy)',{'weight': 'bold', 'size': 15})
g1.set_ylabel('Organ', {'weight': 'bold', 'size': 15})

# stackoverflow post on get_xticklabel being emtpy string: https://stackoverflow.com/questions/11244514/modify-tick-label-text
xtickll = g1.get_xticks().tolist()

g1.set_xticklabels(xtickll, {'weight': 'bold', 'size': 13})
g1.set_yticklabels(g1.get_yticklabels(), {'weight': 'bold', 'size': 13})
plt.tight_layout()
fig_name = '{}/AbsorbedDose_Organs.pdf'.format(data_dir)
plt.savefig(fig_name)
# plt.show()

# g = sns.factorplot(x = 'OrganName_mod', y='ResTime_BqhrPerBq', data=df1, size=6, kind='bar',palette='muted')
# plt.show()
#
# look at relationship btw slope vs RT
df5 = df1[df1['OrganName_mod'] == 'Tumor']

TP_slope_idx = [(1,2),(2,3),(1,3)]
slope_def = ['pIA_{}_{}TP_slope'.format(a,b) for a,b in TP_slope_idx] +\
            ['SUV_{}_{}TP_slope'.format(a,b) for a,b in TP_slope_idx]

for sd in slope_def:
    # get a linear equation
    print(sd)

    # exclue patient 6 and patient 10 when involving time point 3 due to inconsistent time interval
    # the two patients only have 3 time points
    if sd in ['pIA_2_3TP_slope', 'pIA_1_3TP_slope', 'SUV_2_3TP_slope', 'SUV_1_3TP_slope']:
        the_df = df5.ix[~df5['pt_id'].isin(['6','10']),:]
    else:
        the_df = df5

    print(the_df)

    m, b, r_val, p_val, stderr = stats.linregress(the_df['ResTime_BqhrPerBq'], the_df[sd])
    print(m, b, r_val, p_val, stderr)

    sns.set(font='Arial')
    sns.set_context('paper')
    plt.figure(figsize=(11,8.5))

    ax = sns.regplot(x='ResTime_BqhrPerBq', y=sd, data=df5, color='g', line_kws={'label': 'y={0:.3f}x + {1:.3f}, r={2:.3f}'.format(m,b,r_val)}, ci=None)
    ax.legend(fontsize=13)

    # g3 = sns.lmplot(x='ResTime_BqhrPerBq', y='TwoTP_RT_slope', data=df5, ci=None, palette='muted', size=4, scatter_kws={'s':50, 'alpha':1})
    ax.set_xlabel('Residence time (Bq*hr/Bq)', {'weight': 'bold', 'size': 15})
    ax.set_ylabel(sd, {'weight': 'bold', 'size': 15})
    xtickll = ax.get_xticks().tolist()
    ytickll = ax.get_yticks().tolist()
    ytickll_mod = ['{:.3f}'.format(float(ss))for ss in ytickll]

    ax.set_xticklabels(xtickll, {'weight': 'bold', 'size': 13})
    ax.set_yticklabels(ytickll_mod, {'weight': 'bold', 'size': 13})
    plt.tight_layout()
    fig_name = '{}/ResTimeVS{}_Tumors.pdf'.format(data_dir,sd)
    plt.savefig(fig_name, bbox_inches='tight')
    # plt.show()