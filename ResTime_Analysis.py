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
        return 'Adrenal'
    elif s.lower().find('tumor') != -1:
        return 'Tumor'
    else:
        return s

data_dir = '/Users/shuang/Documents/Proj_Neuroblastoma/Data/ResTime_Analysis'

con = mdb.connect('localhost','testuser','test000','UCSFDoseDB')
qr_rt = "SELECT * FROM ResTimeInfo"
df = query_mysql(con, qr_rt)
df['OrganName_mod'] = df['OrganName'].map(srcname_redef)

# save the raw data to csv
csv_fname = '{}/Organ_Tumor_ResTime_2TPslope.csv'.format(data_dir)
df.drop(['id','OrganName_mod'], axis=1).to_csv(csv_fname, index=False)

# selet the df where it's not tumor
df1 = df[(df['OrganName_mod'] != 'Tumor') & (df['OrganName_mod'] != 'TotalBody')]


# TODO: need to fix the plot and make sure it's showing things clearly!
# TODO: PT10 has a total body RT with neg value, .... how do i remedy this? the patient only had 3 imaging time point


sns.set(font='Arial')
sns.set_context('paper')

plt.figure(figsize=(11,8.5))
# f, ax = plt.subplots()
# f.set_size_inches(11,8.5)
g1 = sns.stripplot(x='ResTime_BqhrPerBq', y='OrganName_mod', data=df1, dodge=True, jitter=True, alpha=0.5, zorder=1, palette='dark')
g2 = sns.pointplot(x='ResTime_BqhrPerBq', y='OrganName_mod', data=df1, dodge=0.532, join=False, palette='dark', markers='d', scale=0.75, ci=None)
g1.set_xlim(0,5)
g1.set_xlabel('Residence time (Bq*hr/Bq)',{'weight': 'bold', 'size': 13})
g1.set_ylabel('Organ', {'weight': 'bold', 'size': 13})

# stackoverflow post on get_xticklabel being emtpy string: https://stackoverflow.com/questions/11244514/modify-tick-label-text
xtickll = g1.get_xticks().tolist()

g1.set_xticklabels(xtickll, {'weight': 'bold', 'size': 11})
g1.set_yticklabels(g1.get_yticklabels(), {'weight': 'bold', 'size': 11})
fig_name = '{}/ResTime_Organs.pdf'.format(data_dir)
plt.savefig(fig_name)
# plt.show()

# g = sns.factorplot(x = 'OrganName_mod', y='ResTime_BqhrPerBq', data=df1, size=6, kind='bar',palette='muted')
# plt.show()

# look at relationship btw slope vs RT
df2 = df[df['OrganName_mod'] == 'Tumor']

# get a linear equation
m, b, r_val, p_val, stderr = stats.linregress(df2['ResTime_BqhrPerBq'], df2['TwoTP_RT_slope'])
print(m, b, r_val, p_val, stderr)

sns.set(font='Arial')
sns.set_context('paper')
plt.figure(figsize=(11,8.5))

ax = sns.regplot(x='ResTime_BqhrPerBq', y='TwoTP_RT_slope', data=df2, color='g', line_kws={'label': 'y={0:.3f}x + {1:.3f}, r2={2:.3f}'.format(m,b,r_val**2)}, ci=None)
ax.legend(fontsize=13)

# g3 = sns.lmplot(x='ResTime_BqhrPerBq', y='TwoTP_RT_slope', data=df2, ci=None, palette='muted', size=4, scatter_kws={'s':50, 'alpha':1})
ax.set_xlabel('Residence time (Bq*hr/Bq)', {'weight': 'bold', 'size': 13})
ax.set_ylabel('Two-time-point slope', {'weight': 'bold', 'size': 13})
xtickll = ax.get_xticks().tolist()
ytickll = ax.get_yticks().tolist()
ytickll_mod = ['{:.3f}'.format(float(ss))for ss in ytickll]

ax.set_xticklabels(xtickll, {'weight': 'bold', 'size': 11})
ax.set_yticklabels(ytickll_mod, {'weight': 'bold', 'size': 11})
fig_name = '{}/ResTimeVSTwoTPslope_Tumors.pdf'.format(data_dir)
plt.savefig(fig_name, bbox_inches='tight')
# plt.show()