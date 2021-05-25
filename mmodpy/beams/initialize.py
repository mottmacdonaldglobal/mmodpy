# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:59:31 2020

@author: kevin.stanton
"""

# Import shared libraries
import os
import sys
import pandas as pd
import numpy as np
import shutil
import time

# Import project libraries
from .util import util
from .construct_hysBeams import info_hysBeams
from .iXTRACT_envelope import iXTRACT_envelope
from .Shear_Calcs import Shear_Calcs
from .Autorun import Autorun
from .Mat_Hys_Beam import Mat_Hys_Beam
from .settings import settings
from .construct_elasticBeams import info_elasticBeams
from .Elastic_Beam import Elastic_Beam
from .beamGeometry import beamGeometry

def MHB(iSet=1,
        pInputs=os.path.join(os.getcwd(),'inputs_hysBeams.xlsx'),
        pSave=os.path.join(os.getcwd(),'Output_hysBeams'),
        codeOption='EC2',
        concreteFactor='normal',
        fParams='ASCE 41-17',
        tbOption='default',
        concrete_strain_limits=[
            0.0012,
            0.0024,
            0.0035],
        steel_strain_limits=[
            1.0,
            1.5],
        npoints=401,
        mphisolution='BiSection',
        bilinearflag=False,
        cLimits='explicit',
        cleanCurves=True,
        skipFlag=True,
        pXTRACT=r'C:\Program Files (x86)\TRC\XTRACT\XTRACT.exe',
        mInter=False,
        timingControl='automatic',
        controlMode='default',
        author='user',
        company='Arup',
        job_name='dev',
        showPlots=True,
        write_log=False,
        user_xpj=False,
        settings_path=''):
    
    
    """
    This function intializes an object to enable individual function testing 
    outside of the main pipeline for MAT_HYSTERETIC_BEAM card writing.
    """
    
    
    # Turn off batch mode
    batchMode=False
    
    # Check that inputs file exists and handle if not
    if os.path.isfile(pInputs) == False:
        raise AssertionError(
            'Inputs file not found: ' + str(pInputs))

    # Parse input file
    # Get input file type
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

    # Tag output folder with current time
    timestr = time.strftime("%Y%m%d-%H%M%S")
    pSave = str(pSave + '_' + timestr)
    os.makedirs(pSave)

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
    
    aData = info_hysBeams(data, iRun)
    
    return aData