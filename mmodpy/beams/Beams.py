# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:53:35 2019

@author: kevin.stanton

This script contains functions for the autodyna toolkit.
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


def MAT_HYSTERETIC_BEAM(
        batchMode=False,
        iSet=1,
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
    This function writes MAT_HYSTERETIC_BEAM cards for LS-DYNA. DO NOT TOUCH
    THE MOUSE WHILE CODE IS RUNNING unless the manual intervention option is
    enabled.

    @author: kevin.stanton

    Parameters
    ----------
    batchMode : bool
        Batch mode option
        Default = False -> results only returned for the given iSet
    iSet : int
        Input set for a single analysis (ignored if batchMode = 'on')
    pInputs : str
        Input file path (.txt, .csv, or .xlsx)
        Default = cwd\inputs_hysBeams.xlsx
        Extra values/blank lines are automatically ignored
    pSave : str
        Option to specify a different file path to write results
        Default = cwd\Output_hysBeams
        Specified folder will be tagged with current date/time
    codeOption : str
        'EC2' for Eurocode 2 [default]
        'CSA' for Canadian Standards Association (CSA A23.3-14)
    concreteFactor : str
        Factor to account for lightweight aggregate concrete (or not)
        Only applicable to CSA (known as parameter lambda)
        'normal' (lamda=1) [default]; 'semi-light' (lambda=0.85); or 'light'
        (lambda=0.75)
    fParams : str or bool
        Options to compute/write plastic rotation thresholds
        Valid arguments = 'ASCE 41-17' (default), 'ASCE 41-13', or False
        If False, no element deletion will take place
        PRSi/PRTi = 0 if required input data is missing
    tbOption : str
        Treatment of T-Beams ('default' or 'in-slab')
        'in-slab' -> density_reduced = density*(Area_eff/Area_total)
    concrete_strain_limits : list of floats
        Limiting concrete strains for enveloping PM curve
        For more refined analyses use:
            [0.0012,0.0015,0.0018,0.0021,0.0024,0.0027,0.003,0.0033,0.0035]
        ASCE 41-17, Section 10.3.3.1 [Usable Strain Limits]:
            [0.002,0.003]
    steel_strain_limits : list of floats
        Limiting steel strains for enveloping PM curve
        Defined as: (yield stress / Young's modulus) x steel_strain_limits[i]
        default = [1.0, 1.5, 2.0]
    npoints : int
        Number of points desired on the final, enveloped PM curve
        Default = 401
    mphisolution : str
        Solution method for moment curvature plot
        default = 'BiSection'
        'MUDN' -> iteration at Minimum Unbalanced Displacement Norm
    bilinearflag : bool
        Flag to use bi-linearized moment curvature data in LS-DYNA
        Default = False
    cLimits : str
        Yield moment interpretation ('strict', 'flex', 'explicit')
        'strict' and 'flex' are ML methods
        'explicit' (default) is based on checking developed strains against
        user-defined limits
    cleanCurves:
        Flag to write smoothed moment curvature data to LS-DYNA
        default = True
    skipFlag : True
        Flag to skip XTRACT analysis if appropriate data already exists
        Relevant data must be saved in working directory ([ID].txt files)
        Default = True
    pXTRACT : str
        XTRACT path
        default = r'C:\Program Files (x86)\TRC\XTRACT\XTRACT.exe'
        The code automatically checks the os for XTRACT but sometimes
        environment settings force the path to be defined explicitely
    mIter : bool
        Option to enable manual intervention mode for XTRACT analysis
        Default = False
    timingControl : str
        Toggle for automatic (CPU-based) analysis timing or fixed timing.
        Acceptable arguments:
            'automatic' -> CPU usage based analysis timing [default]
            'fixed' -> fixed analysis time from 'settings.py'
    controlMode : str
        Mouse movement control options ('default' or 'mod1')
        If 'default' causes an error, try 'mod1'
    author : str
        Name of author
        Default = 'user'
    company : str
        Company name
        Default = 'Arup'
    job_name : str
        Project name
        Default = 'dev'
    showPlots : bool
        Show XTRACT plots when running section analysis (True or False)
        May run faster if 'False'
        Default = True
    write_log : bool
        Option to write the contents of the terminal to a log file
        Default = False
    user_xpj : bool
        Option to use a user-defined .xpj to run XTRACT
        File name must match info_filename for given iSet and file must be 
        stored in the cwd
    settings_path : str
        Path to settings file for user-defined mouse click positions
        Default setting used if unspecified

    Returns
    -------
    Beam-RC_Beam_[ID].key : file(s)
        LS-DYNA keyword file is written for each beam instance
        Includes section and material card data
    autobeams.key : file
        LS-DYNA keyword containing all batch results in one file
    [ID].txt : file
        A text file containing the raw XTRACT output is saved for each
        beam instance
    M_phi_xx_[ID].png, M_phi_yy_[ID].png : files
        Images of the interpretted moment-curvature plots for both loading
        directions are saved each beam instance
    Moment_Curvature_[ID].csv : file
        Processed moment curvature data for beam instance [ID]
    Moment_Curvature_RAW_[ID].csv : file
        Raw moment curvature data for beam instance [ID]
    Moment_Interaction_[ID].csv : file
        Axial load-moment interaction data for beam instance [ID]
        Applies to last entries in steel and concrete strain limit lists
    Moment_Interaction[ID].csv : file
        Axial load-moment interaction data for beam instance [ID]
        Includes all results for the different steel and concrete strain limits
    Moment_Interaction_Envelope[ID].csv : file
        Enveloped axial load-moment interaction data for beam instance [ID]
    Moment_Plastic_Rotation_[ID].csv : file
        Moment-plastic rotation interaction data for beam instance [ID]
    PMxx[ID].csv, PMyy[ID].csv : file
        Enveloped axial load-moment interaction data for beam instance [ID]
    jsonData : str
        Contents of Beam-RC_Beam_[ID].key written in JSON format
    keyPath : str
        Path to generated file 'autobeams.key'
    hysBeams.log : file (optional)
        The console output written to a log file
        Not written in write_log = False
    """

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

    # Add pSave to path
    sys.path.insert(0, pSave)
    
    # Store output in a log file if needed
    if write_log:
        class Logger(object):
            def __init__(self):
                self.terminal = sys.stdout
                self.log = open(os.path.join(pSave,'hysBeams.log'), 'w')

            def write(self, message):
                self.terminal.write(message)
                self.log.write(message)  

            def flush(self):
                self.log.flush() 

            def close(self):
                self.log.close()                
        sys.stdout = Logger()

    # Run required analyses
    if batchMode:

        keyNames = []
        iRun = 0
        fRun = iData.shape[1]  # Total number of runs for print statements
        while iRun < iData.shape[1]:
            # Initialize analysis inputs
            aData = info_hysBeams(data, iRun)

            # Process mouse control settings
            mcSettings = settings(aData)

            # Print warning statement if 'BH' is not included in the ID in
            # batchData.csv
            util.getName(aData)

            if envelopeflag:

                iSet = iRun + 1
                print('\n--------------------')
                print('Beam ' + aData.file_name +
                      ' (' + str(iSet) + ' of ' + str(fRun) + ')')
                print('--------------------\n')
                print('   Section type: ' + str(aData.section_type) + '\n')

                if XTRACT == 'on':
                    if not user_xpj:
                        print('   Creating XTRACT input file... \n')
                        saveName = str(aData.file_name) + '.txt'
                        if skipFlag:
                            if os.path.isfile(saveName) == False:
                                # Create the section in XTRACT
                                section_area = iXTRACT_envelope(aData).write_EnvelopeXtract(
                                    aData, concrete_strain_limits, steel_strain_limits, mphisolution)
                                # Run XTRACT
                                Autorun(
                                    aData,
                                    mcSettings).autorun_xtract(
                                    aData,
                                    n_strainlim=len(concrete_strain_limits) *
                                    len(steel_strain_limits))
                            else:
                                shutil.copy(saveName, pSave)
                        else:
                            # Move results to current output folder
                            shutil.copy(saveName, pSave + '\\' +
                                        str(aData.file_name) + '.txt')
                            # Create the section in XTRACT
                            section_area = iXTRACT_envelope(aData).write_EnvelopeXtract(
                                aData, concrete_strain_limits, steel_strain_limits, mphisolution)
                            # Run XTRACT
                            Autorun(
                                aData,
                                mcSettings).autorun_xtract(
                                aData,
                                n_strainlim=len(concrete_strain_limits) *
                                len(steel_strain_limits))
                    else:
                        xpj_name = str(aData.file_name) + '.xpj'
                        if os.path.isfile(xpj_name) == False:
                            raise AssertionError('Cannot find user-defined .xpj file.')
                        else:
                            print('   User-defined .xpj file located. Attempting to analyze in XTRACT... ')
                            # Run XTRACT
                            Autorun(
                                aData,
                                mcSettings, xpj_name).autorun_xtract(
                                aData,
                                n_strainlim=len(concrete_strain_limits) *
                                len(steel_strain_limits))

                # Calculate the PM interaction envelope
                util.getPM(
                    aData,
                    len(concrete_strain_limits) *
                    len(steel_strain_limits),
                    npoints)

                # Shear calculations
                print('   Performing shear calculations... \n')
                shear = Shear_Calcs(aData)
                Ned = float(aData.N_load)
                Vrdx, Vrdy = shear.Shear_cap(Ned)
                PRS, PRT = shear.prLimits()

                # Create individual .key files
                print('   Writing keyword data...')
                card = Mat_Hys_Beam(
                    aData,
                    aData.section_area,
                    envelopeflag,
                    len(concrete_strain_limits) *
                    len(steel_strain_limits),
                    bilinearflag)
                json_data = card.write_mat(aData, PRS, PRT)

                # Populate list of .key file names
                keyNames.append(os.path.join(aData.pSave,
                                             util.getName(aData, '.key')))
			
			# Update counter
            iRun = iRun + 1

        # Combine .key data into one file (autoBeams.key)
        util.combineKeys(keyNames)
        # Move autoBeams.key to output directory
        shutil.move('autobeams.key', aData.pSave)

        # Report on completion
        print('-- Finished --')
        
        # Close logger
        if write_log:
            Logger().close()
            
        # Return keyword path if pSave specified, JSON data if not
        if pSave:
            keyPath = os.path.join(aData.pSave,'autobeams.key')
            return keyPath
        else:
            return json_data

    elif not batchMode:

        # Specify data set to run on
        iRun = iSet - 1

        # Initialize analysis inputs
        aData = info_hysBeams(data, iRun)

        # Process mouse control settings
        mcSettings = settings(aData)

        if envelopeflag:

            print('\n--------------------')
            print('Beam ' + aData.file_name)
            print('--------------------\n')
            print('   Section type: ' + str(aData.section_type) + '\n')

            if XTRACT == 'on':
                if not user_xpj:
                    print('   Creating XTRACT input file... \n')
                    saveName = os.getcwd() + '\\' + str(aData.file_name) + '.txt'
                    if skipFlag:
                        if os.path.isfile(saveName) == False:
                            # Create the section in XTRACT
                            section_area = iXTRACT_envelope(aData).write_EnvelopeXtract(
                                aData, concrete_strain_limits, steel_strain_limits, mphisolution)
                            # Run XTRACT
                            Autorun(
                                aData,
                                mcSettings).autorun_xtract(
                                aData,
                                n_strainlim=len(concrete_strain_limits) *
                                len(steel_strain_limits))
                        else:
                            shutil.copy(saveName, pSave)
                    else:
                        # Move results to current output folder
                        shutil.copy(saveName, pSave + '\\' +
                                    str(aData.file_name) + '.txt')
                        # Create the section in XTRACT
                        section_area = iXTRACT_envelope(aData).write_EnvelopeXtract(
                            aData, concrete_strain_limits, steel_strain_limits, mphisolution)
                        # Run XTRACT
                        Autorun(
                            aData,
                            mcSettings).autorun_xtract(
                            aData,
                            n_strainlim=len(concrete_strain_limits) *
                            len(steel_strain_limits))
                else:
                        xpj_name = str(aData.file_name) + '.xpj'
                        if os.path.isfile(xpj_name) == False:
                            raise AssertionError('Cannot find user-defined .xpj file.')
                        else:
                            print('   User-defined .xpj file located. Attempting to analyze in XTRACT... ')
                            # Run XTRACT
                            Autorun(
                                aData,
                                mcSettings, xpj_name).autorun_xtract(
                                aData,
                                n_strainlim=len(concrete_strain_limits) *
                                len(steel_strain_limits))

            # Calculate PM interaction envelope
            util.getPM(aData,
                       len(concrete_strain_limits) *
                       len(steel_strain_limits),
                       npoints)

            # Calculate shear capacity
            print('   Performing shear calculations... \n')
            shear = Shear_Calcs(aData)
            Ned = float(aData.N_load)
            Vrdx, Vrdy = shear.Shear_cap(Ned)
            PRS, PRT = shear.prLimits()

            # Create .key file
            print('   Writing keyword data...')
            card = Mat_Hys_Beam(
                aData,
                aData.section_area,
                envelopeflag,
                len(concrete_strain_limits) *
                len(steel_strain_limits),
                bilinearflag)
            json_data = card.write_mat(aData, PRS, PRT)

        # Report on completion
        print('-- Finished --')

        # Close logger
        if write_log:
            Logger().close()
            
        # Return keyword path if pSave specified, JSON data if not
        if pSave:
            keyPath = os.path.join(aData.pSave, util.getName(aData, '.key'))
            return keyPath
        else:
            return json_data


def MAT_ELASTIC(
        batchMode=False,
        iSet=1,
        pInputs=os.path.join(os.getcwd(), 'inputs_elasticBeams.xlsx'),
        pSave=os.path.join(os.getcwd(), 'Output_elasticBeams'),
        tbOption='default',
        author='user',
        company='Arup',
        job_name='dev',
        write_log=False):
    """
    This function writes MAT_ELASTIC and section cards for LS-DYNA beam
    elements.

    @author: kevin.stanton

    Parameters
    ----------
    batchMode : bool
        Batch mode option ('on' or 'off')
        Default = False -> results only returned for the given iSet
    iSet : int
        Input set for a single analysis (ignored if batchMode = 'on')
    pInputs : str
        Input file path (.txt, .csv, or .xlsx)
        Default = os.getcwd() + '\\inputs.xlsx'
        Extra values/blank lines are automatically ignored
        May alternatively be defined as a list of strings in lieu of all other
        arguments
    pSave : str
        Option to specify a different file path to write results
        Default = os.getcwd() + '\\Output_elasticBeams'
        Specified folder will be tagged with current date/time
    tbOption : str
        Treatment of T-Beams ('default' or 'in-slab')
        'in-slab' -> density_reduced = density*(Area_eff/Area_total)
    author : str
        Name of author
        Default = 'user'
    company : str
        Company name
        Default = 'Arup'
    job_name : str
        Project name
        Default = 'dev'
    write_log : bool
        Option to write the contents of the terminal to a log file
        Default = False

    Returns
    -------
    Beam-RC_Beam_[ID].key : file(s)
        An LS-DYNA keyword file is written to the output directory for each
        beam analyzed
        Includes section and material card data
    autobeams.key : file
        LS-DYNA keyword containing all batch results in one file is written to
        the output directory
    jsonData : str
        Contents of Beam-RC_Beam_[ID].key written in JSON format
    keyPath : str
        Path to generated file 'autobeams.key'
    elasticBeams.log : file (optional)
        The console output is written to a log file the pSave directory stamped
        with the analysis start date/time
        Not written if write_log = False
    """

    # Specify interpreter order
    iData_order = ['info_filename',
                   'Ec',
                   'den',
                   'v',
                   'section_type',
                   'cover',
                   'diameter_circ',
                   'rect_width',
                   'rect_height',
                   't_width',
                   't_height',
                   'f_width',
                   'f_thick',
                   'bX',
                   'bY',
                   'tX',
                   'tY',
                   'cX',
                   'cY']

    # Parse inputs and handle the case where pInputs is a list of strings
    if type(pInputs) == list:
        batchMode = False
        iSet = 1
        pSave = False
        tbOption = 'default'
        author = 'user'
        company = 'Arup'
        job_name = 'dev'
        iData = [iData_order, pInputs]
        iData = pd.DataFrame(np.transpose(np.array(iData)))
    else:
        # Check that inputs file exists and handle if not
        if os.path.isfile(pInputs) == False:
            raise AssertionError(
                'Inputs file not found: ' + str(pInputs))
        # Get input file type
        ext = os.path.splitext(pInputs)[1]
        # Read data
        if ext == '.txt':
            iData = pd.read_csv(
                pInputs,
                dtype=str,
                header=None,
                sep=r',\s+',
                engine='python')
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

    if not pSave == False:
        # Tag output folder with current time
        timestr = time.strftime("%Y%m%d-%H%M%S")
        pSave = pSave + '_' + timestr
        os.makedirs(pSave)
        # Add pSave to path
        sys.path.insert(0, pSave)
    
    # Store output in a log file if needed
    if write_log:
        class Logger(object):
            def __init__(self):
                self.terminal = sys.stdout
                self.log = open(os.path.join(pSave,'elasticBeams.log'), 'w')

            def write(self, message):
                self.terminal.write(message)
                self.log.write(message)  

            def flush(self):
                self.log.flush()  
                
            def close(self):
                self.log.close()
        sys.stdout = Logger()        

    # Append static inputs to data array
    sData = [
        author,
        company,
        job_name,
        batchMode,
        iSet,
        tbOption,
        pSave,
        write_log]
    data = iData, sData

    # Run required analyses
    if batchMode:

        keyNames = []
        iRun = 0
        while iRun < iData.shape[1]:

            # Initialize analysis inputs
            aData = info_elasticBeams(data, iRun)

            # Print warning statement if 'BH' is not included in the ID in
            # batchData.csv
            util.getName(aData)

            print('\n--------------')
            print('Beam ' + aData.file_name)
            print('-------------- \n')
            print('   Section type: ' + str(aData.section_type) + '\n')

            # Compute geometric section data
            section_area = beamGeometry().runElastic(aData)

            # Create individual .key files
            print('   Writing keyword data...')
            card = Elastic_Beam(aData)
            json_data = card.write_mat(aData)

            # Populate list of .key file names
            keyNames.append(os.path.join(aData.pSave,
                                         util.getName(aData, '.key')))

            # Update counter
            iRun = iRun + 1

        # Combine .key data into one file (autoBeams.key)
        util.combineKeys(keyNames)
        # Move autoBeams.key to output directory
        shutil.move('autobeams.key', aData.pSave)

        # Close logger
        if write_log:
            Logger().close()
            
        # Return keyword path if pSave specified, JSON data if not
        if pSave:
            keyPath = os.path.join(aData.pSave,'autobeams.key')
            return keyPath
        else:
            return json_data

    elif not batchMode:

        # Specify data set to run on
        iRun = iSet - 1

        # Initialize analysis inputs
        aData = info_elasticBeams(data, iRun)

        print('\n--------------')
        print('Beam ' + aData.file_name)
        print('-------------- \n')
        print('   Section type: ' + str(aData.section_type) + '\n')

        # Compute geometric section data
        section_area = beamGeometry().runElastic(aData)

        # Create .key file
        print('   Writing keyword data...')
        card = Elastic_Beam(aData)
        json_data = card.write_mat(aData)
        
        # Close logger
        if write_log:
            Logger().close()
            
        # Return keyword path if pSave specified, JSON data if not
        if pSave:
            keyPath = os.path.join(aData.pSave, util.getName(aData, '.key'))
            return keyPath
        else:
            return json_data
            
