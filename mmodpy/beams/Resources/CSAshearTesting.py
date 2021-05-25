# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:58:30 2020

@author: kevin.stanton
"""


import autodyna
import os
import sys
import pandas as pd
import numpy as np
import shutil
import time
from autodyna.construct_hysBeams import info_hysBeams
from autodyna.Shear_Calcs_test import Shear_Calcs

# Initialize
pInputs='C:\\Users\\kevin.stanton\\Documents\\inputs_hysBeams_028.xlsx'
settings_path = 'C:\\ArupLocal\\GitLab\\AutomatedCard_Generator_LS-DYNA\\autodyna\\settings.txt'
batchMode=False
iSet=45
pInputs=pInputs
codeOption='CSA'
write_log=True
concreteFactor='normal'
fParams='ASCE 41-17'
tbOption='default'
concrete_strain_limits=[
    0.0012,
    0.0024,
    0.0035]
steel_strain_limits=[
    1.0,
    1.5]
npoints=401
mphisolution='BiSection'
bilinearflag=False
cLimits='explicit'
cleanCurves=True
skipFlag=True
pXTRACT=r'C:\Program Files (x86)\TRC\XTRACT\XTRACT.exe'
mInter=False
timingControl='automatic'
controlMode='default'
author='user'
company='Arup'
job_name='dev'
showPlots=True
user_xpj=False
pSave=os.path.join(os.getcwd(),'Output_hysBeams')
ext = os.path.splitext(pInputs)[1]
# Specify interpreter order
iData_order = ['info_filename',
               'Fc',
               'Ft',
               'ey',
               'ecu',
               'esp',
               'ef',
               'Ec',
               'den',
               'v',
               'Fy',
               'Fu',
               'esh',
               'esu',
               'Es',
               'section_type',
               'cover',
               'diameter_circ',
               'bar_diam_circ',
               'num_bars_circ',
               'rect_width',
               'rect_height',
               'n_t',
               'n_b',
               'db_t',
               'db_b',
               'num_bars_col',
               'bar_diam_col',
               't_width',
               't_height',
               'f_width',
               'f_thick',
               'f_lay',
               'f_bar',
               'f_spac',
               't_nt',
               't_tbar',
               't_nb',
               't_bbar',
               'bX',
               'bY',
               'tX',
               'tY',
               'cX',
               'cY',
               'cBar_centerD',
               'cBar_vFlangeD',
               'cBar_hFlangeD',
               'cBar_vFlangeN',
               'cBar_hFlangeN',
               'cBar_rows',
               'M_load',
               'N_load',
               'V_load',
               'd_shear',
               's_shear',
               'n_legs',
               'gc_mat',
               'gs_mat',
               'Fywd',
               'phl_xx',
               'phl_yy',
               'rb_condition',
               'c_condition',
               'ss_ratio_s',
               'ss_ratio_t']

# Read data
if ext == '.txt':
    iData = pd.read_csv(pInputs, dtype=str, header=None,
                        sep=',\s+', engine='python')
elif ext == '.csv':
    iData = pd.read_csv(pInputs, dtype=str, header=None)
elif ext == '.xlsx':
    iData = pd.read_excel(pInputs, dtype=str, header=None)
# Remove whitespaces from strings
iData = iData.applymap(lambda x: x.strip() if isinstance(x, str) else x)
# Remove data rows not spacified in iData_order
iData = iData[iData.iloc[:, 0].isin(iData_order)]
# Remove second column if full input format detected
if 'List' in iData.iloc[0, 1]:
    iData = iData.drop(columns=1)
# Re-index
iData = iData.set_index(0)
# Sort according to iData_order
iData = iData.loc[iData_order]
# Re-index columns to match iRun handles
iData.columns = pd.RangeIndex(0, len(iData.columns))

# Append static inputs to data array
XTRACT = 'on'
devMode = 'on'
envelopeflag = True  # use PM enveloping approach
sData = [
    author,
    company,
    job_name,
    pXTRACT,
    mInter,
    batchMode,
    iSet,
    fParams,
    tbOption,
    devMode,
    cLimits,
    cleanCurves,
    showPlots,
    controlMode,
    codeOption,
    concreteFactor,
    timingControl,
    pSave,
    write_log,
    settings_path]
data = iData, sData
iRun = iSet - 1

###########################################
############# Important Stuff #############
###########################################

#### Initialize analysis inputs ####
aData = info_hysBeams(data, iRun)

#### Initialize shear capacity class ####
shear = Shear_Calcs(aData)

#### Define explicit inputs for shear calculations
Ned = float(aData.N_load)


#### Calculate shear capacity in the x direction (x = t = width) ####
width = 254
height = 939.8
n_s = 2
d_s = 12.7
s_spacing = 304.8
n_t = 3
d_t = 22.225
n_b = 2
d_b = 22.225
Vrdx = shear.Shear_cap(Ned, width, height, n_s, d_s, s_spacing, n_t, d_t, n_b, d_b)
PRS = shear.prLimits()

#### Calculate shear capacity in the x direction (x = t = width) ####
height = 254
width = 939.8
n_s = 2
d_s = 12.7
s_spacing = 304.8
n_t = 2
d_t = 22.225
n_b = 3
d_b = 22.225
Vrdy = shear.Shear_cap(Ned, width, height, n_s, d_s, s_spacing, n_t, d_t, n_b, d_b)
PRT = shear.prLimits()