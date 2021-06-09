# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 14:22:33 2020

@author: kevin.stanton

This project is a WIP.
"""

import os
import pandas as pd
import numpy as np

def run_TCF(inp_path, thf_path, this_path):
    '''
    This function executes an Oasys T-HIS FastTCF script on the local machine.

    Parameters
    ----------     
    inp_path : str
        Full file path to the desired .inp file 
    thf_path : str
        Full file path to the .thf file containing the data to be processed
    this_path : str
        Full file path to the Oasys T-HIS executable on the local machine       

    Returns
    -------
    files : files(s)
        List of new or modified files produced in the thf_path directory
   '''
    
    # Import libraries
    import subprocess
    import time
    import os
    
    # Initialize
    cwd = os.getcwd()
    thf_dir = os.path.dirname(thf_path)
    os.chdir(thf_dir)
    start_time = time.time()
    
    # Define commands for automatically running FAST-TCF   
    com1 = '"' + os.path.splitext(str(this_path))[0] + '"'
    com2 = ' ' + str(thf_path) + ' '
    com3 = " -tcf=" + '"' + str(inp_path) + '"'
    com4 = " -batch -noconsole -exit"
    cmds = com1+com2+com3+com4
    # Create string to write to batch file
    batchData = 'net session>nul 2>&1\n'+\
    '\n'+\
    ':main\n'+\
    str(cmds)+'\n'+\
    'exit'
    # Write batch data and copy to appropriate directory
    with open('thfBatch.bat','w') as file:
        file.write(batchData)
    # Export model results
    p = subprocess.Popen('thfBatch.bat')
    p.communicate()
    
    # Record new/modified file paths
    files = []
    for dirname,subdirs,fnames in os.walk(thf_dir):
        for fname in fnames:
            full_path = os.path.join(dirname, fname)
            mtime = os.stat(full_path).st_mtime
            if mtime > start_time:
                files.append(fname)
    
    # Revert to original cwd
    os.chdir(cwd)
    
    return files
        
def find_replace(walk_dir, find, replace, filePattern):
    '''
    This function walk through a directory and replaces specific text within 
    files matching the filePattern extension.
    
    Parameters
    ----------     
    walk_dir : str
        Directory within which to find and replace text
    find : str
        Text string to be replaced
    replace : str
        Text string to replace with   
    filePattern : str
        File extension(s) whithin which to search/replace (e.g. "*.txt")

    Returns
    -------
    Modified file(s) : files(s)
        All changes are made by overwriting existing files
    '''
    
    import os
    import fnmatch
    
    walk_dir = os.path.normpath(walk_dir)
    for path, dirs, files in os.walk(os.path.abspath(walk_dir)):
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            with open(filepath) as f:
                s = f.read()
            s = s.replace(find, replace)
            with open(filepath, "w") as f:
                f.write(s)

def run_dyna(exe_path, key_path, max_cpu=100, overwrite = False, lic_path=False):
    '''
    This function executes an LS-DYNA analysis on the local machine.

    Parameters
    ----------     
    exe_path : str
        Full file path to the desired LS-DYNA executable   
    key_path : str
        Full file path to the LS-DYNA keyword file to be run   
    max_cpu : int
        Maximum number of CPUs
    overwrite : bool (optional)
        Option to overwrite results if they exit (no effect if not)   
    lic_path : bool or str (optional)
        Path to local license file (ignore if using network license)

    Returns
    -------
    files : files(s)
        List of new or modified files produced in the thf_path directory
   '''

    # Import libraries
    import os
    import subprocess
    import multiprocessing
    import time
    
    # Define variables
    cwd = os.getcwd()
    exe_path = os.path.normpath(exe_path)
    key_path = os.path.normpath(key_path)
    key_dir = os.path.dirname(key_path)
    os.chdir(key_dir)
    nthreads = min(str(min(int(multiprocessing.cpu_count()/2),4)),max_cpu)
    start_time = time.time()
    
    # Run LS-DYNA analysis
    if overwrite:
        owFlag = 1
    elif not os.path.isfile(key_dir + '\\messag'):
        owFlag = 1
    elif not 'N o r m a l    t e r m i n a t i o n' in open(key_dir + '\\messag').read():
        owFlag = 1
    else:
        owFlag = 0
    if owFlag == 1:
        com1 = '"' + os.path.splitext(str(exe_path))[0] + '" '
        com2 = 'i='+ key_dir + '.key '
        com3 = 'g='+ key_dir + '.ptf '
        com4 = 'f='+ key_dir + '.thf '
        com5 = 'o='+ key_dir + '.otf '
        com6 = 'ncpu='+nthreads 
        if lic_path:
            com7 = 'set LSTC_LICENSE=local\nset LSTC_FILE=' + lic_path + '\n'
            cmds = com1+com2+com3+com4+com5+com6+com7
        else:
            cmds = com1+com2+com3+com4+com5+com6
        # Create string to write to batch file        
        batchData = 'net session>nul 2>&1\n'+\
        '\n'+\
        ':main\n'+\
        str(cmds)+'\n'+\
        'exit'
        # Write batch data
        with open('dynaBatch.bat','w') as file:
            file.write(batchData)
        # Run LS-DYNA automatically for all motions
        p = subprocess.Popen('dynaBatch.bat')
        p.communicate()
 
    # Record new/modified file paths
    files = []
    for dirname,subdirs,fnames in os.walk(key_dir):
        for fname in fnames:
            full_path = os.path.join(dirname, fname)
            mtime = os.stat(full_path).st_mtime
            if mtime > start_time:
                files.append(fname)
    
    # Revert to original cwd
    os.chdir(cwd)
    
    return files
        
def lin_scale(seeds, targets, tMinMax=(0, 10), weights=pd.DataFrame()):
    """
    Linearly scales seed motions to match a specified target response
    spectra. The scale factor is determined according to a weighting scheme either
    specified by the user or assumed to be constant for all periods.

    =========
    Arguments
    =========
    
    seeds : DataFrame
        - Multidimensional DataFrame of seed periods and spectral ordinates
        - Must contain columns titles denoting the name of the source motion,
        component ( or ), and ordinate/abscissa flag (_T or _Sa)

    targets : DataFrame
        - Target response spectra (x = period, y = spectral acceleration)

    tMinMax : tuple [optional, default=(0, 10)]
        - If specified, equal weights are applied to periods within Tmin and Tmax.
        - This input is overridden if user-defined weight are specified.

    weights : 2D array or dict [optional]
        - Discrete weighting factors (x = period, y = weighting factor). The sum of
        all specified weighting factors must equal 1. Seed/target periods outside
        the specified range will not be scaled. This input takes precedence over
        tMinMax if it is specified.
        
    =======
    Returns
    =======
    
        regData_raw : DataFrame
            - Concatenated seeds and targets regularized and unscaled   
            
        regData_scaled : DataFrame
            - Concatenated seeds and targets regularized and scaled 
            
        sfOptimized : Numpy Array
            - Optimized scale factors for all periods
    """
    
    # Import libraries
    import pandas as pd
    import numpy as np
    from scipy.optimize import minimize
    
    # Disable new pandas warnings
    pd.options.mode.chained_assignment = None

    #############################
    # Parse and regularize data #
    #############################

    # Format seeds and targets
    seeds = pd.DataFrame(seeds)
    targets = pd.DataFrame(targets)
    seeds.columns = seeds.iloc[0]
    targets.columns = targets.iloc[0]
    seeds = seeds.reindex(seeds.index.drop(0))
    targets = targets.reindex(targets.index.drop(0))
    seeds = seeds.dropna()
    targets = targets.dropna()
    seeds = seeds.reindex(sorted(seeds.columns), axis=1)
    targets = targets.reindex(sorted(targets.columns), axis=1)
    seeds = seeds.apply(pd.to_numeric, errors='ignore')
    targets = targets.apply(pd.to_numeric, errors='ignore')
    tMinMax = (tMinMax[0], tMinMax[1])
    weights = pd.DataFrame(weights)

    # Identify governing min/max periods from seed and target data
    cData = pd.concat([targets, seeds], axis=1)
    rawPeriods = cData.loc[:, cData.columns.str.contains('_T')]
    minT = min(rawPeriods.min())
    maxT = max(rawPeriods.max())

    # Crop specified periods of interest outside the seed/target range
    tMin = max(min(tMinMax), minT)
    tMax = min(max(tMinMax), maxT)

    # Define regularized periods
    regLen = max(50, len(seeds), len(targets))
    T = np.linspace(minT, maxT, regLen)
    T = [round(float(T[i]), 3) for i in range(0, len(T))]

    # Create weighting array from tMinMax if not specified by the user
    if weights.empty:
        weights = pd.DataFrame({'T': [0, max(0, tMin - 0.001), tMin, tMax, min(
            10, tMax + 0.001), max(10, tMax + 0.002)], 'wf': [0, 0, 0.5, 0.5, 0, 0]})
    else:
        # Format user-defined weighting array
        weights.columns = ['T', 'wf']
        weight_Tmin = min(weights['T'])
        weight_Tmax = max(weights['T'])
        wBoundsLow = []
        wBoundsHigh = []
        wBoundsLow.insert(0, {'T': max(0, weight_Tmin - 0.0001), 'wf': 0})
        wBoundsLow.insert(0, {'T': 0, 'wf': 0})
        wBoundsHigh.insert(0, {'T': max(10, weight_Tmax + 0.002), 'wf': 0})
        wBoundsHigh.insert(0, {'T': min(10, weight_Tmax + 0.001), 'wf': 0})
        weights = pd.concat([pd.DataFrame(wBoundsLow), weights,
                             pd.DataFrame(wBoundsHigh)], ignore_index=True)

    # Interpolate spectral accelerations and weights for regularized periods
    regData = pd.DataFrame({'T': T})
    wScale = sum(np.interp(T, weights['T'], weights['wf']))
    regData['weights'] = np.interp(T, weights['T'], weights['wf'] / wScale)

    for col in cData.columns:
        cInd = cData.columns.get_loc(col)
        if cData.columns.str.contains('_Sa')[cInd]:
            regData[col] = np.interp(
                T, cData[cData.columns[cInd + 1]], cData[col])
    
    # Drop residual nan values from regData
    regData = regData.dropna()
    
    #########
    # Scale #
    #########

    # Function to compute room mean squared error for the prescribed weighting
    def RMSE(sf, *args):

        # Compute weighted mean squared residuals
        seed1, target1, weight1 = args
        wMSR = target1.copy()
        for i in range(0, len(target)):
            wMSR[i] = weight1[i] * (target1[i] - seed1[i] * sf)**2

        # Compute RMSE
        RMSE = sum(wMSR)

        return RMSE

    # Minimize RMSE between targets and all seeds
    regData_raw = regData.copy()
    regData_scaled = regData.copy()  # Place holder
    regData_copy = regData.copy()  # Prevent chained modifications
    sfOptimized = np.arange(len(regData.columns) - 3)
    sfOptimized = [float(sfOptimized[i]) for i in sfOptimized]
    targetI = regData_copy[regData_copy.columns[2]]
    weight = regData_copy[regData_copy.columns[1]]
    sfMin = 0.1
    sfMax = 5
    for i in range(3, len(regData.columns)):

        # Compute the error vector
        seed = regData_copy[regData_copy.columns[i]]
        target = targetI.copy()
        bnds = tuple([sfMin, sfMax])
        bnds = tuple([bnds])
        res = minimize(
            RMSE,
            x0=(sfMin),
            args=(
                seed,
                target,
                weight),
            bounds=bnds)
        sfOptimized[i - 3] = res.x

        # Compute the weighted seeds
        regData_scaled[regData_scaled.columns[i]
                       ] = regData[regData.columns[i]] * res.x

    return regData_raw, regData_scaled, sfOptimized

def write_sra_column(input_xlsx, th_dir, out_dir = os.getcwd(), viscd=0.001):
    '''
    --WIP--
    
    This function creates the main .key file (MODEL.key) that is used to run
    SRA in LS-DYNA for all motions in the input directory (located one level
    below mwd). Returns a node ID at the surface for post-processing. 
    '''
    
    import math
    import init_geo
    import time

    # Initialize objects
    params = [input_xlsx, th_dir]
    data = init_geo.update(params)
    
    # Check that input time histories exist
    if not os.path.isdir(th_dir):
        raise AssertionError('Directory for input time-histories not found.')
    
    # Create analysis folders
    dir_name = 'sra_' + str(int(time.time()))
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    os.makedirs(os.path.join(out_dir,dir_name))
    os.makedirs(os.path.join(out_dir,dir_name,'DYNA'))
    os.makedirs(os.path.join(out_dir,dir_name,'DYNA','Time_Histories_(vel)'))
    os.makedirs(os.path.join(out_dir,dir_name,'DYNA','Images'))
        
    # Initialize variables
    if data.back_gen == 'Generated':
        strn, gg0, shst, damp = _read_back(data.exn)
    else:
        strn, gg0, shst, damp = _backbone(data.exn)
    s_id = data.s_id
    viscd = data.vis_damp
    bed_sh = data.bed_sh
    bed_cp = data.bed_cp
    bed_de = data.bed_de
    
    # Create output folder
    loc = out_dir + '\\DYNA\\SRA'
    if os.path.isdir(loc) == False:
        os.makedirs(loc)

    # Lysmer properties
    vs_bed = bed_sh  # Bedrock shear wave velocity (m/sec)
    vp_bed = bed_cp  # Bedrock shear wave velocity (m/sec)
    dens_bot = bed_de  # Bedrock density (kg/m3)
    vp = bed_cp  # Define p-wave velocity (m/sec)

    # Define label ids
    s_id_10 = str(s_id).rjust(10)
    lys_id = str(s_id + len(data)).rjust(10)
    lys_id8 = str(s_id + len(data)).rjust(8)
    nhis_id = str(s_id)
    nhis_10 = str(s_id).rjust(10)
    spc_id = str(s_id + 1).rjust(10)

    # Define control cards
    lctm_id = str(s_id + len(data) + 1).rjust(10)

    # Write the data into the key file
    loc = out_dir + '\\DYNA\\SRA\\sra.key'
    fwrite = open(loc, 'w')
    fwrite.write('*KEYWORD\n')
    fwrite.write('$\n')

    # Define control cards
    fwrite.write('$ =============\n')
    fwrite.write('$ CONTROL cards\n')
    fwrite.write('$ =============\n')
    fwrite.write('$\n')
    fwrite.write('*CONTROL_ENERGY\n')
    fwrite.write('$$    HGEN      RWEN    SLNTEN     RYLEN\n')
    fwrite.write('         2         0         2         2\n')
    fwrite.write('*CONTROL_SOLID\n')
    fwrite.write('         1         0         0         0         0' +
                 '         0\n')
    fwrite.write('*CONTROL_SOLUTION\n')
    fwrite.write('         0         0         0     15000\n')
    fwrite.write('*CONTROL_TIMESTEP\n')
    fwrite.write('$$  DTINIT    TSSFAC      ISDO    TSLIMT     DT2MS\
      LCTM     ERODE     MS1ST\n')
    fwrite.write('       0.0       0.8         0       0.0       0.0'
                 + lctm_id + '         0         0\n')
    fwrite.write('$\n')

    # Define curve for timestep
    z_10 = '         0'
    fwrite.write('*DEFINE_CURVE_TITLE\n')
    fwrite.write('Max timestep\n')
    fwrite.write(lctm_id + z_10 + z_10 + z_10 + z_10 + z_10 + z_10 + '\n')
    fwrite.write('                 0.0              2.0E-4\n')
    fwrite.write('               100.0              2.0E-4\n')
    fwrite.write('$\n')

    # Define database cards
    fwrite.write('$ ==============\n')
    fwrite.write('$ DATABASE cards\n')
    fwrite.write('$ ==============\n')
    fwrite.write('$\n')
    fwrite.write('*DATABASE_BINARY_D3PLOT\n')
    fwrite.write('    1.0E-1         0         0         0         0\n')
    fwrite.write('*DATABASE_BINARY_D3THDT\n')
    fwrite.write('    5.0E-3         0\n')
    fwrite.write('$\n')
    fwrite.write('*DATABASE_EXTENT_BINARY\n')
    fwrite.write('$$   NEIPH     NEIPS    MAXINT    STRFLG    SIGFLG    EPSFLG\
        RLTFLG    ENGFLG\n')
    fwrite.write('         0         0         0         1         0'
                 + '         0         0         0\n')
    fwrite.write('$$  CMPFLG    IEVERP    BEAMIP     DCOMP      SHGE     STSSZ\
        N3THDT    IALEMAT\n')
    fwrite.write('         0         0         0         0         0'
                 + '         0         0         0\n')
    fwrite.write('$$ NINTSLD   PKP_SEN      SCLP     HYDRO     MSSCL     THERM\
        INTOUT    NODOUT\n')
    fwrite.write('         0         0       0.0         0         0'
                 + '         0                    \n')
    fwrite.write('$$    DTDT\n')
    fwrite.write('         0         0         0\n')
    fwrite.write('*DATABASE_FORMAT\n')
    fwrite.write('         0         1\n')
    fwrite.write('$\n')
    fwrite.write('$\n')

    # Define hourglass
    fwrite.write('$ ===============\n')
    fwrite.write('$ HOURGLASS cards\n')
    fwrite.write('$ ===============\n')
    fwrite.write('$\n')
    fwrite.write('*HOURGLASS\n')
    fwrite.write(s_id_10 + '         4    1.0E-3         0'
                 + '       0.0       0.0       0.0       0.0\n')
    fwrite.write('$\n')

    # Write the section cards
    fwrite.write('$ =============\n')
    fwrite.write('$ SECTION cards\n')
    fwrite.write('$ =============\n')
    fwrite.write('$\n')
    fwrite.write('*SECTION_SOLID_TITLE\n')
    fwrite.write('Soil\n')
    fwrite.write(s_id_10 + '         1         0\n')
    fwrite.write('$\n')
    fwrite.write('*SECTION_DISCRETE_TITLE\n')
    fwrite.write('Lysmer Damper\n')
    fwrite.write(lys_id + '         0       0.0       0.0       0.0'
                 + '       0.0\n')
    fwrite.write('       0.0       0.0\n')
    fwrite.write('$\n')

    # Write the material cards
    mat_ids = range(s_id, s_id + len(data))
    mat_ids = [str(item).rjust(10) for item in mat_ids]
    dens = (data.inp['gamma (kN/m3)'] / 0.00981).values
    bulk_mod = vp**2 * dens - 4.0 / 3.0 * data.inp['G0 (kPa)'].values
    dens = np.around(dens, decimals=2)
    dens = [str('%.2E' % item).rjust(10) for item in dens]
    bulk_mod = np.around(bulk_mod, decimals=2)
    bulk_mod = [str('%.2E' % item) for item in bulk_mod]
    bulk_mod = [item.replace("+", "") for item in bulk_mod]
    bulk_mod = [str(item).rjust(10) for item in bulk_mod]
    cut_off = -1.0E10 * np.ones(len(mat_ids))
    cut_off = [str('%.2E' % item) for item in cut_off]
    cut_off = [item.replace("+", "") for item in cut_off]
    cut_off = [str(item).rjust(10) for item in cut_off]
    lc_ids = np.arange(s_id, s_id + len(mat_ids))
    lc_ids = [str(item).rjust(10) for item in lc_ids]
    ones_10 = len(mat_ids) * ["       1.0"]
    zeros_10 = len(mat_ids) * ["         0"]
    null_10 = len(mat_ids) * ["          "]
    mat_wr1 = [list(i) for i in zip(mat_ids, dens, bulk_mod, cut_off,
                                    zeros_10, ones_10, zeros_10, zeros_10)]
    mat_wr1 = [''.join(mat_wr1[i]) for i in range(len(mat_wr1))]
    mat_wr2 = [list(i) for i in zip(zeros_10, ones_10, lc_ids, ones_10,
                                    zeros_10, zeros_10, zeros_10, zeros_10)]
    mat_wr2 = [''.join(mat_wr2[i]) for i in range(len(mat_wr2))]
    mat_wr3 = [list(i) for i in zip(zeros_10, zeros_10, zeros_10, zeros_10,
                                    zeros_10, null_10, zeros_10, zeros_10)]
    mat_wr3 = [''.join(mat_wr3[i]) for i in range(len(mat_wr3))]
    mat_wr4 = [list(i) for i in zip(zeros_10, zeros_10, zeros_10, zeros_10,
                                    zeros_10, null_10, zeros_10, zeros_10)]
    mat_wr4 = [''.join(mat_wr4[i]) for i in range(len(mat_wr4))]
    fwrite.write('$ ====================\n')
    fwrite.write('$ MAT (Material) cards\n')
    fwrite.write('$ ====================\n')
    fwrite.write('$\n')
    fwrite.write('*MAT_DAMPER_VISCOUS\n')
    fwrite.write(lys_id + '       1.0\n')
    fwrite.write('$\n')
    for i in range(len(mat_ids)):
        fwrite.write('*MAT_HYSTERETIC_SOIL_TITLE\n')
        fwrite.write('Layer ' + str(i + 1) + ' - ' +
                     str(data.inp['Soil Type'][i]) + '\n')
        fwrite.write(mat_wr1[i])
        fwrite.write('\n')
        fwrite.write(mat_wr2[i])
        fwrite.write('\n')
        fwrite.write(mat_wr3[i])
        fwrite.write('\n')
        fwrite.write(mat_wr4[i])
        fwrite.write('\n')
        fwrite.write('$\n')

    # Write strength curve cards
    curve_id = range(s_id, s_id + len(data))
    c_tag = [str(item).rjust(10) for item in curve_id]
    zer_10, sfa, sfo = '         0', '     0.866', '     0.866'
    fwrite.write('$\n')
    for i in range(len(data)):
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Layer ' + str(i + 1) + ' - ' +
                     str(data.inp['Soil Type'][i]) + '\n')
        fwrite.write(
            c_tag[i] +
            zer_10 +
            sfa +
            sfo +
            zer_10 +
            zer_10 +
            zer_10 +
            '\n')
        str_lvl_temp = [str("%.2E" % (item / 100.0)).rjust(20) for item in
                        strn.tolist()]
        sh_str_temp = [str("%.6E" % item).rjust(20)
                       for item in shst[i].tolist()]
        for j in range(len(strn)):
            fwrite.write(str_lvl_temp[j] + sh_str_temp[j] + '\n')
        fwrite.write('$\n')

    # Write the part cards
    part_ids = range(s_id, s_id + len(data))
    part_ids = [[i] for i in part_ids]
    part_ids = [item for sublist in part_ids for item in sublist]
    part_ids = [str(item).rjust(10) for item in part_ids]
    part_name = np.arange(1, len(part_ids) + 1)
    part_name = ['Layer ' +
                 str(part_name[i]) +
                 ' - ' +
                 str(data.inp['Soil Type'][i]) for i in range(len(part_name))]
    mat_ids = range(s_id, s_id + len(data))
    mat_ids = [str(item).rjust(10) for item in mat_ids]
    sec_id = s_id_10
    hg_id = s_id_10
    zeroes_10 = "         0"
    fwrite.write('$ ==========\n')
    fwrite.write('$ PART cards\n')
    fwrite.write('$ ==========\n')
    fwrite.write('$\n')
    for i in range(len(part_ids)):
        fwrite.write('*PART\n')
        fwrite.write(part_name[i])
        fwrite.write('\n')
        t_str = str(data.inp['Soil Type'][i])
        if t_str == "Clay" or t_str == "Silty Clay":
            fwrite.write('$PR_PART_COL ' + str(part_ids[i]) + ' 7c1ee0f0')
        elif t_str == "Peat":
            fwrite.write('$PR_PART_COL ' + str(part_ids[i]) + ' 792fe0f0')
        elif t_str == "Silt":
            fwrite.write('$PR_PART_COL ' + str(part_ids[i]) + ' 7957e0f0')
        elif t_str == "Sand":
            fwrite.write('$PR_PART_COL ' + str(part_ids[i]) + ' 7bd960f0')
        else:
            fwrite.write('$PR_PART_COL ' + str(part_ids[i]) + ' 7f5960f0')
        fwrite.write('\n')
        fwrite.write(part_ids[i] + sec_id + mat_ids[i] + zeroes_10 + hg_id +
                     zeroes_10 + zeroes_10 + zeroes_10 + "\n")
    fwrite.write('$\n')
    fwrite.write('*PART\n')
    fwrite.write('$PR_PART_COL ' + str(str(s_id + len(data))) + ' 7dbbe0f0')
    fwrite.write('\n')
    fwrite.write('Lysmer dampers\n')
    fwrite.write(lys_id + lys_id + lys_id + '         0         0\
             0         0         0\n')
    fwrite.write('$\n')

    # Write the node cards
    node_ids = range(s_id, s_id + (2 * 2 * (len(data) + 1) + 1) + 3)
    node_ids = [str(item).rjust(8) for item in node_ids]
    z_elev = data.inp['Top Depth (m)'].values
    z_elev = np.append(z_elev, float(data.inp['Bottom Depth (m)'].iloc[-1]))
    z_elev = np.around(-z_elev, decimals=2)
    z_elev = z_elev.tolist()
    z_elev.append(min(z_elev) - 1)
    z_elev = [str(item).rjust(16) for item in z_elev]
    col_wid = 1
    x_coord = [0, col_wid]
    y_coord = [0, col_wid]
    x_coord = [str(item).rjust(16) for item in x_coord]
    y_coord = [str(item).rjust(16) for item in y_coord]
    node_ids = [[i] for i in node_ids]
    node_x = 2 * (len(data) + 2) * [[i] for i in x_coord]
    node_y = (len(data) + 2) * [2 * [i] for i in y_coord]
    node_y = [item for sublist in node_y for item in sublist]
    node_y = [[i] for i in node_y]
    node_z = [2 * 2 * [i] for i in z_elev]
    node_z = [item for sublist in node_z for item in sublist]
    node_z = [[i] for i in node_z]
    zer_8 = len(node_ids) * ["        0"]
    nl = len(node_ids) * ["\n"]
    node_wr = [list(i) for i in zip(node_ids, node_x, node_y,
                                    node_z, zer_8, zer_8, nl)]

    def flatten(l): return flatten(l[0]) + (flatten(l[1:])
                                            if len(l) > 1 else [])\
        if isinstance(l, list) else [l]
    node_wr = [flatten(node_wr[i]) for i in range(len(node_wr))]
    node_wr = [''.join(node_wr[i]) for i in range(len(node_wr))]
    fwrite.write('$ ==========\n')
    fwrite.write('$ NODE cards\n')
    fwrite.write('$ ==========\n')
    fwrite.write('$\n')
    fwrite.write('*NODE\n')
    [fwrite.write(node_wr[i]) for i in range(len(node_wr))]
    fwrite.write('$\n')

    # Write the element cards
    solid_ids = range(s_id, s_id + (1 * 1 * (len(data)) + 1))
    solid_ids = [str(item).rjust(8) for item in solid_ids]
    part_ids = range(s_id, s_id + ((len(data))))
    part_ids = [1 * 1 * [i] for i in part_ids]
    part_ids = [item for sublist in part_ids for item in sublist]
    part_ids = [str(item).rjust(8) for item in part_ids]
    node_ids = np.array(range(s_id, s_id + 2 * 2 * (len(data) + 1)))
    max_node_id = node_ids[-1]
    z_elev = data.inp['Top Depth (m)'].values
    z_elev = np.append(z_elev, float(data.inp['Bottom Depth (m)'].iloc[-1]))
    z_elev = np.around(-z_elev, decimals=2)
    z_elev = z_elev.tolist()
    z_elev = [str(item).rjust(16) for item in z_elev]
    node_ids = node_ids.reshape((len(z_elev), 2, 2))
    node_1 = np.delete(node_ids, 0, axis=1)
    node_1 = np.delete(node_1, -1, axis=2)
    node_1 = np.delete(node_1, -1, axis=0)
    node_1 = node_1.reshape((-1, 1))
    node_1 = [item for sublist in node_1 for item in sublist]
    node_1 = [str(item).rjust(8) for item in node_1]
    node_2 = np.delete(node_ids, 0, axis=1)
    node_2 = np.delete(node_2, 0, axis=2)
    node_2 = np.delete(node_2, -1, axis=0)
    node_2 = node_2.reshape((-1, 1))
    node_2 = [item for sublist in node_2 for item in sublist]
    node_2 = [str(item).rjust(8) for item in node_2]
    node_3 = np.delete(node_ids, -1, axis=1)
    node_3 = np.delete(node_3, 0, axis=2)
    node_3 = np.delete(node_3, -1, axis=0)
    node_3 = node_3.reshape((-1, 1))
    node_3 = [item for sublist in node_3 for item in sublist]
    node_3 = [str(item).rjust(8) for item in node_3]
    node_4 = np.delete(node_ids, -1, axis=1)
    node_4 = np.delete(node_4, -1, axis=2)
    node_4 = np.delete(node_4, -1, axis=0)
    node_4 = node_4.reshape((-1, 1))
    node_4 = [item for sublist in node_4 for item in sublist]
    node_4 = [str(item).rjust(8) for item in node_4]
    node_5 = np.delete(node_ids, 0, axis=1)
    node_5 = np.delete(node_5, -1, axis=2)
    node_5 = np.delete(node_5, 0, axis=0)
    node_5 = node_5.reshape((-1, 1))
    node_5 = [item for sublist in node_5 for item in sublist]
    node_5 = [str(item).rjust(8) for item in node_5]
    node_6 = np.delete(node_ids, 0, axis=1)
    node_6 = np.delete(node_6, 0, axis=2)
    node_6 = np.delete(node_6, 0, axis=0)
    node_6 = node_6.reshape((-1, 1))
    node_6 = [item for sublist in node_6 for item in sublist]
    node_6 = [str(item).rjust(8) for item in node_6]
    node_7 = np.delete(node_ids, -1, axis=1)
    node_7 = np.delete(node_7, 0, axis=2)
    node_7 = np.delete(node_7, 0, axis=0)
    node_7 = node_7.reshape((-1, 1))
    node_7 = [item for sublist in node_7 for item in sublist]
    node_7 = [str(item).rjust(8) for item in node_7]
    node_8 = np.delete(node_ids, -1, axis=1)
    node_8 = np.delete(node_8, -1, axis=2)
    node_8 = np.delete(node_8, 0, axis=0)
    node_8 = node_8.reshape((-1, 1)).tolist()
    node_8 = [item for sublist in node_8 for item in sublist]
    node_8 = [str(item).rjust(8) for item in node_8]
    nl = len(solid_ids) * ["\n"]
    solid_wr = [
        list(i) for i in zip(
            solid_ids,
            part_ids,
            node_6,
            node_5,
            node_8,
            node_7,
            node_2,
            node_1,
            node_4,
            node_3,
            nl)]
    solid_wr = [flatten(solid_wr[i]) for i in range(len(solid_wr))]
    solid_wr = [''.join(solid_wr[i]) for i in range(len(solid_wr))]
    nod_dam = range(s_id, s_id + (2 * 2 * (len(data) + 1) + 1) + 3)
    ndu1, ndu2, ndu3, ndu4 = nod_dam[-8], nod_dam[-7], nod_dam[-6], nod_dam[-5]
    ndd1, ndd2, ndd3, ndd4 = nod_dam[-4], nod_dam[-3], nod_dam[-2], nod_dam[-1]
    ndu1, ndu2 = str(ndu1).rjust(8), str(ndu2).rjust(8)
    ndu3, ndu4 = str(ndu3).rjust(8), str(ndu4).rjust(8)
    ndd1, ndd2 = str(ndd1).rjust(8), str(ndd2).rjust(8)
    ndd3, ndd4 = str(ndd3).rjust(8), str(ndd4).rjust(8)
    dens = (data.inp['gamma (kN/m3)'] / 0.00981).values
    dens = [float(str('%.2E' % item)) for item in dens]
    damp_hor = str(
        float(col_wid) /
        2 *
        float(col_wid) /
        2 *
        vs_bed *
        dens_bot).rjust(16)
    damp_ver = str(
        float(col_wid) /
        2 *
        float(col_wid) /
        2 *
        vp_bed *
        dens_bot).rjust(16)
    sd_idx = str(s_id + 0).rjust(8)
    sd_idy = str(s_id + 1).rjust(8)
    sd_idz = str(s_id + 2).rjust(8)
    fwrite.write('$ =============\n')
    fwrite.write('$ ELEMENT cards\n')
    fwrite.write('$ =============\n')
    fwrite.write('$\n')
    fwrite.write('*ELEMENT_SOLID\n')
    [fwrite.write(solid_wr[i]) for i in range(len(solid_wr))]
    fwrite.write('$\n')
    fwrite.write('*ELEMENT_DISCRETE\n')
    fwrite.write(str(s_id + 0).rjust(8) + lys_id8 + ndu1 +
                 ndd1 + sd_idx + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 1).rjust(8) + lys_id8 + ndu1 +
                 ndd1 + sd_idy + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 2).rjust(8) + lys_id8 + ndu1 +
                 ndd1 + sd_idz + damp_ver + '       0     0.0\n')
    fwrite.write(str(s_id + 3).rjust(8) + lys_id8 + ndu2 +
                 ndd2 + sd_idx + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 4).rjust(8) + lys_id8 + ndu2 +
                 ndd2 + sd_idy + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 5).rjust(8) + lys_id8 + ndu2 +
                 ndd2 + sd_idz + damp_ver + '       0     0.0\n')
    fwrite.write(str(s_id + 6).rjust(8) + lys_id8 + ndu3 +
                 ndd3 + sd_idx + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 7).rjust(8) + lys_id8 + ndu3 +
                 ndd3 + sd_idy + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 8).rjust(8) + lys_id8 + ndu3 +
                 ndd3 + sd_idz + damp_ver + '       0     0.0\n')
    fwrite.write(str(s_id + 9).rjust(8) + lys_id8 + ndu4 +
                 ndd4 + sd_idx + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 10).rjust(8) + lys_id8 + ndu4 +
                 ndd4 + sd_idy + damp_hor + '       0     0.0\n')
    fwrite.write(str(s_id + 11).rjust(8) + lys_id8 + ndu4 +
                 ndd4 + sd_idz + damp_ver + '       0     0.0\n')
    fwrite.write('$\n')

    # Define define cards
    fwrite.write('$ ============\n')
    fwrite.write('$ DEFINE cards\n')
    fwrite.write('$ ============\n')
    fwrite.write('$\n')
    fwrite.write('*DEFINE_SD_ORIENTATION\n')
    fwrite.write(str(s_id + 0).rjust(10) + '         0       1.0       0.0'
                 + '       0.0\n')
    fwrite.write(str(s_id + 1).rjust(10) + '         0       0.0       1.0'
                 + '       0.0\n')
    fwrite.write(str(s_id + 2).rjust(10) + '         0       0.0       0.0'
                 + '       1.0\n')
    fwrite.write('$\n')

    # Define boundary SPC cards
    ndd1, ndd2, ndd3, ndd4 = nod_dam[-4], nod_dam[-3], nod_dam[-2], nod_dam[-1]
    ndd1, ndd2 = str(ndd1).rjust(10), str(ndd2).rjust(10)
    ndd3, ndd4 = str(ndd3).rjust(10), str(ndd4).rjust(10)
    fwrite.write('$ ==============\n')
    fwrite.write('$ BOUNDARY cards\n')
    fwrite.write('$ ==============\n')
    fwrite.write('$\n')
    fwrite.write('*BOUNDARY_SPC_SET\n')
    fwrite.write(spc_id + '         0         1         1         1'
                 + '         0         0         0\n')
    fwrite.write('$\n')
    fwrite.write('$: Cross-reference summary for SET_NODE ' + s_id_10 + '\n')
    fwrite.write('$:------------------------------------------\n')
    fwrite.write('$: SPC <No label>\n')
    fwrite.write('$\n')
    fwrite.write('*SET_NODE_LIST_TITLE\n')
    fwrite.write('Nodes on bottom plane\n')
    fwrite.write(spc_id + '       0.0       0.0       0.0       0.0\n')
    fwrite.write(ndd1 + ndd2 + ndd3 + ndd4 + '\n')
    fwrite.write('$\n')

    # Define the constrained node sets
    node_ids = range(s_id, s_id + 2 * 2 * (len(data) + 1))
    z_elev = data.inp['Top Depth (m)'].values
    z_elev = np.append(z_elev, float(data.inp['Bottom Depth (m)'].iloc[-1]))
    z_elev = np.around(-z_elev, decimals=2)
    z_elev = z_elev.tolist()
    col_wid = 1
    x_coord = [0, col_wid]
    y_coord = [0, col_wid]
    node_x = 2 * (len(data) + 1) * [[i] for i in x_coord]
    node_x = [float(node_x[i][0]) for i in range(len(node_x))]
    node_y = (len(data) + 1) * [2 * [i] for i in y_coord]
    node_y = [item for sublist in node_y for item in sublist]
    node_z = [2 * 2 * [i] for i in z_elev]
    node_z = [item for sublist in node_z for item in sublist]
    node_ids = np.array(node_ids).reshape(-1, 1)
    node_x = np.array(node_x).reshape(-1, 1)
    node_y = np.array(node_y).reshape(-1, 1)
    node_z = np.array(node_z).reshape(-1, 1)
    node_com = np.concatenate((node_ids, node_x, node_y, node_z), axis=1)
    per_nodes = []
    for i in range(len(z_elev)):
        z_log = node_com[:, 3] == z_elev[i]
        x1_log = node_com[:, 1] == 1
        x2_log = node_com[:, 1] == 0
        y1_log = node_com[:, 2] == 1
        y2_log = node_com[:, 2] == 0
        log = z_log * (x1_log + x2_log + y1_log + y2_log)
        per_nodes.append(node_com[log][:, 0].tolist())
    for i in range(len(per_nodes)):
        per_nodes[i] = map(int, per_nodes[i])
        per_nodes[i] = map(str, per_nodes[i])
        per_nodes[i] = [str(item).rjust(10) for item in per_nodes[i]]
        per_nodes[i] = [per_nodes[i][x:x + 8]
                        for x in range(0, len(per_nodes[i]), 8)]
    for i in range(len(per_nodes)):
        per_nodes[i] = [''.join(per_nodes[i][j])
                        for j in range(len(per_nodes[i]))]
    base_nodes = node_com[node_com[:, 3] == z_elev[-1]][:, 0]
    base_nodes = list(map(int, base_nodes))
    base_nodes = list(map(str, base_nodes))
    base_nodes = [str(item).rjust(10) for item in base_nodes]
    base_nodes = [base_nodes[x:x + 8] for x in range(0, len(base_nodes), 8)]
    base_nodes = [''.join(base_nodes[j]) for j in range(len(base_nodes))]
    set_ids = np.arange(
        s_id + len(data) + 1,
        s_id + len(z_elev) + len(data) + 2)
    set_ids = [str(item).rjust(10) for item in set_ids]
    set_ids_s = np.arange(
        s_id +
        len(data) +
        1,
        s_id +
        len(z_elev) +
        len(data) +
        2).tolist()
    set_ids_s = list(map(str, set_ids_s))
    ndu1, ndu2, ndu3, ndu4 = nod_dam[-8], nod_dam[-7], nod_dam[-6], nod_dam[-5]
    ndu1, ndu2 = str(ndu1).rjust(10), str(ndu2).rjust(10)
    ndu3, ndu4 = str(ndu3).rjust(10), str(ndu4).rjust(10)
    zer = "         0"
    fwrite.write('$\n')
    fwrite.write('$ =================\n')
    fwrite.write('$ CONSTRAINED cards\n')
    fwrite.write('$ =================\n')
    fwrite.write('$\n')
    for i in range(len(z_elev)):
        fwrite.write('*CONSTRAINED_NODAL_RIGID_BODY_SPC\n')
        fwrite.write(
            set_ids[i] +
            zer +
            set_ids[i] +
            zer +
            zer +
            zer +
            zer +
            '\n')
        fwrite.write('       1.0       0.0       7.0\n')
        fwrite.write('$\n')
        fwrite.write('$: Cross-reference summary for SET_NODE ' +
                     set_ids_s[i] + '\n')
        fwrite.write('$:---------------------------------------\n')
        fwrite.write('$: NODAL_RIGID_BODY <No label>\n')
        fwrite.write('$\n')
        fwrite.write('*SET_NODE_LIST_TITLE\n')
        fwrite.write('Side nodes at z=' + str(z_elev[i]) + '\n')
        fwrite.write(set_ids[i] + zer + zer + zer + zer + '\n')
        [fwrite.write(per_nodes[i][j] + '\n')
         for j in range(len(per_nodes[i]))]
        fwrite.write('$\n')

    # Define the damping cards
    vs_av = np.average(data.inp['Vs (m/s)'].values,
                       weights=data.inp['Thickness (m)'].values)
    height = np.sum(data.inp['Thickness (m)'].values)
    omeg_n = 1 / (4 * height / vs_av) * 2 * math.pi
    beta = round(max(-viscd * 2 / omeg_n, -0.02 / (math.pi * 20.0)), 6)
    beta = str("%.3E" % beta).rjust(10)
    fwrite.write('$\n')
    fwrite.write('$ =============\n')
    fwrite.write('$ DAMPING cards\n')
    fwrite.write('$ =============\n')
    fwrite.write('$\n')
    fwrite.write('*DAMPING_PART_STIFFNESS_SET\n')
    fwrite.write(s_id_10 + beta + '\n')
    fwrite.write('$\n')

    # Defind the set cards
    part_ids = range(s_id, s_id + ((len(data))))
    part_ids = [str(item).rjust(10) for item in part_ids]
    part_wr_in, part_wr_tot, counter = [], [], 0
    while counter < len(part_ids):
        for i in range(8):
            try:
                part_wr_in.append(part_ids[counter + i])
            except BaseException:
                continue
        part_wr_in.append('\n')
        counter = counter + 8
        part_wr_tot.append(part_wr_in)
        part_wr_in = []
    part_wr_tot = [''.join(part_wr_tot[i]) for i in range(len(part_wr_tot))]
    fwrite.write('$ =========\n')
    fwrite.write('$ SET cards\n')
    fwrite.write('$ =========\n')
    fwrite.write('$\n')
    fwrite.write('$: Cross-reference summary for SET_PART 1\n')
    fwrite.write('$:---------------------------------------\n')
    fwrite.write('$: DAMPING <No label>\n')
    fwrite.write('$\n')
    fwrite.write('*SET_PART_LIST_TITLE\n')
    fwrite.write('Soil Set Part\n')
    fwrite.write(s_id_10 + '       0.0       0.0       0.0       0.0\n')
    [fwrite.write(part_wr_tot[i])for i in range(len(part_wr_tot))]
    fwrite.write('$\n')

    # Define the load cards
    ndu1, ndu2, ndu3, ndu4 = nod_dam[-8], nod_dam[-7], nod_dam[-6], nod_dam[-5]
    ndu1, ndu2 = str(ndu1).rjust(10), str(ndu2).rjust(10)
    ndu3, ndu4 = str(ndu3).rjust(10), str(ndu4).rjust(10)
    damp_hor = str(
        float(col_wid) /
        2 *
        float(col_wid) /
        2 *
        vs_bed *
        dens_bot).rjust(10)
    damp_ver = str(
        float(col_wid) /
        2 *
        float(col_wid) /
        2 *
        1500 *
        dens_bot).rjust(10)
    fwrite.write('$\n')
    fwrite.write('$ ==========\n')
    fwrite.write('$ LOAD cards\n')
    fwrite.write('$ ==========\n')
    fwrite.write('$\n')
    fwrite.write('*LOAD_NODE_POINT\n')
    fwrite.write(ndu1 + '         1        11' + damp_hor + '         0\n')
    fwrite.write(ndu1 + '         2        12' + damp_hor + '         0\n')
    fwrite.write(ndu1 + '         3        13' + damp_ver + '         0\n')
    fwrite.write(ndu2 + '         1        11' + damp_hor + '         0\n')
    fwrite.write(ndu2 + '         2        12' + damp_hor + '         0\n')
    fwrite.write(ndu2 + '         3        13' + damp_ver + '         0\n')
    fwrite.write(ndu3 + '         1        11' + damp_hor + '         0\n')
    fwrite.write(ndu3 + '         2        12' + damp_hor + '         0\n')
    fwrite.write(ndu3 + '         3        13' + damp_ver + '         0\n')
    fwrite.write(ndu4 + '         1        11' + damp_hor + '         0\n')
    fwrite.write(ndu4 + '         2        12' + damp_hor + '         0\n')
    fwrite.write(ndu4 + '         3        13' + damp_ver + '         0\n')
    fwrite.write('$\n')

    # Write the database history nodes
    node_ids = range(s_id, s_id + 2 * 2 * (len(data) + 1))
    node_ids = [str(item).rjust(10) for item in node_ids]
    node_db = np.array(node_ids).reshape(-1,
                                         1)[(node_x == 0) & (node_y == 0)].tolist()
    fwrite.write('$: Cross-reference summary for SET_NODE ' + nhis_id + '\n')
    fwrite.write('$:------------------------------------------\n')
    fwrite.write('$: DATABASE_HISTORY <No label>\n')
    fwrite.write('$\n')
    fwrite.write('*SET_NODE_LIST\n')
    fwrite.write(nhis_10 + '       0.0       0.0       0.0       0.0\n')
    nod_wr_in, nod_wr_tot, counter = [], [], 0
    while counter < len(node_db):
        for i in range(8):
            try:
                nod_wr_in.append(node_db[counter + i])
            except BaseException:
                continue
        nod_wr_in.append('\n')
        counter = counter + 8
        nod_wr_tot.append(nod_wr_in)
        nod_wr_in = []
    nod_wr_tot = [''.join(nod_wr_tot[i]) for i in range(len(nod_wr_tot))]
    [fwrite.write(nod_wr_tot[i])for i in range(len(nod_wr_tot))]
    fwrite.write('\n')
    fwrite.write('$\n')
    fwrite.write('*DATABASE_HISTORY_NODE_SET\n')
    fwrite.write(nhis_10 + '\n')
    fwrite.write('$\n')

    # Write the database history solids
    solid_ids = range(s_id, s_id + len(data))
    solid_ids = [str(item).rjust(10) for item in solid_ids]
    solid_db = solid_ids
    fwrite.write('$: Cross-reference summary for SET_SOLID' +
                 str(s_id) + '\n')
    fwrite.write('$:------------------------------------------\n')
    fwrite.write('$: DATABASE_HISTORY <No label>\n')
    fwrite.write('$\n')
    fwrite.write('*SET_SOLID\n')
    fwrite.write(s_id_10 + '\n')
    solid_wr_in, solid_wr_tot, counter = [], [], 0
    while counter < len(solid_db):
        for i in range(8):
            try:
                solid_wr_in.append(solid_db[counter + i])
            except BaseException:
                continue
        solid_wr_in.append('\n')
        counter = counter + 8
        solid_wr_tot.append(solid_wr_in)
        solid_wr_in = []
    solid_wr_tot = [''.join(solid_wr_tot[i])
                    for i in range(len(solid_wr_tot))]
    [fwrite.write(solid_wr_tot[i])for i in range(len(solid_wr_tot))]
    fwrite.write('\n')
    fwrite.write('$\n')
    fwrite.write('*DATABASE_HISTORY_SOLID_SET\n')
    fwrite.write(s_id_10 + '\n')
    fwrite.write('$\n')

    # Close the write handle
    fwrite.write('*END\n')
    fwrite.close()

    # Return output data
    return max_node_id

def th_to_key(input_dir, out_dir, filePattern='*.csv'):
    '''
    This function writes input velocity time histories keyword for LS-DYNA
    using acceleration time-history data (m/sec2) from a csv file.
    
    Parameters
    ----------
        input_dir : str
            Directory containing acceleration time history data {m/sec2}
            Will only read .csv files formatted for the program SIREN
            If results are to be used with SirenAutoRun.exe, file titles must
            be formatted as TH[id].csv, where [id] is odd for x-dir and even 
            for y-dir          
        out_dir : str
            Directory for storing LS-DYNA keyword data           
        filePattern : str
            Extension for files containing time history data (.csv is the 
            only current option)

    Returns
    -------    
        th_accel : DataFrame
            Collection of all acceleration time histories found in input_dir             
        th_vel : DataFrame
            Collection of all velocity time histories computed from th_accel            
        [gm_data_files].key : file(s)
            Keywords files with velocity time history load curve data for two 
            loading directions per file
    '''
    
    from math import Decimal
    import shutil
    import fnmatch
    from scipy import integrate
                           
    # Clear output directory if it already exists and create it if not
    if not os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    # Populate acceleration time history data
    th_accel = pd.DataFrame()
    for path, dirs, files in os.walk(os.path.abspath(input_dir)):
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            th_accel_temp = pd.read_csv(filepath)
            name = list(th_accel_temp.columns.values)[0]
            names = [name + ': time', name + ': accel']
            th_accel_temp.columns = names
            th_accel_temp = th_accel_temp[~th_accel_temp.isin(
                ['Time', 'Acceleration', 'TIME', 'ACCEL'])]
            th_accel_temp = th_accel_temp.dropna()
            th_accel_temp = th_accel_temp.reset_index(drop=True)
            th_accel = pd.concat([th_accel, th_accel_temp], axis=1)
            th_accel = th_accel.astype(float)

    # Convert accelerations to velocity
    th_vel = th_accel
    for i in range(0, len(th_accel.columns)):
        if (i % 2) == 0:
            th_vel.iloc[:, i] = th_accel.iloc[:, i]
        else:
            th_vel.iloc[:, i] = integrate.cumtrapz(
                th_accel.iloc[:, i], th_accel.iloc[:, i - 1], initial=0)

    # Write velocty data to keywords
    for i in range(0, int(len(th_vel.columns) / 2 - 1)):
        termination = str(round(max(th_vel.iloc[:, i * 2]), 1)).rjust(10)
        rec_name = 'REC' + '_' + str(i + 1) + 'x' + str(i + 2) + 'y'
        fName = out_dir + '\\REC' + '_' + \
            str(i + 1) + 'x' + str(i + 2) + 'y' + '_VEL_1G.key'
        fwrite = open(fName, 'w')
        fwrite.write('*KEYWORD\n')
        fwrite.write('$\n')
        fwrite.write('*INCLUDE\n')
        fwrite.write('..\MODEL.key\n')
        fwrite.write('$\n')
        fwrite.write('*CONTROL_TERMINATION\n')
        fwrite.write(
            termination +
            '         0       0.0       0.0       0.0         0\n')
        fwrite.write('$\n')
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write(rec_name + '_X_VEL\n')
        fwrite.write(
            '        11         0       1.0       1.0       0.0       0.0         0         0\n')
        for j in range(len(th_vel.iloc[:, i * 2])):
            fwrite.write('%s%s\n' %
                         (str('%.4E' %
                              Decimal(str(th_vel.iloc[j, i * 2]))).rjust(20), str('%.4E' %
                                                                                  Decimal(str(th_vel.iloc[j, i * 2 + 1]))).rjust(20)))
        fwrite.write('$\n')
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write(rec_name + '_Y_VEL\n')
        fwrite.write(
            '        12         0       1.0       1.0       0.0       0.0         0         0\n')
        for j in range(len(th_vel.iloc[:, i * 2 + 2])):
            fwrite.write('%s%s\n' %
                         (str('%.4E' %
                              Decimal(str(th_vel.iloc[j, i * 2 + 2]))).rjust(20), str('%.4E' %
                                                                                      Decimal(str(th_vel.iloc[j, i * 2 + 3]))).rjust(20)))
        fwrite.write('$\n')
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write(rec_name + '_Z_VEL\n')
        fwrite.write(
            '        13         0       1.0       1.0       0.0       0.0         0         0\n')
        fwrite.write('      0.00000000E+00      0.00000000E+00\n')
        fwrite.write('      5.00000000E+02      0.00000000E+00\n')
        fwrite.write('$\n')
        fwrite.write('*END')
        fwrite.close()

    return th_accel, th_vel

def _read_back(input_file_path):

    '''
    This is a private function that reads soil backbone curves.

    Parameters
    ----------
    input_file_path : str or path
        Full file path to input data spreadsheet 
        See example on GitHub for formatting

    Returns
    -------
    main.key : file
        LS-DYNA keyword data files
    '''
    
    import pandas as pd
    
    # Read soil degradation curves and profile data from the input file
    back = pd.read_excel(input_file_path, sheet_name='Backbone')
    inp = pd.read_excel(input_file_path, sheet_name='Soil', header=0)

    # Read the data for the defined strain levels
    strn = np.array(back.iloc[0:, 0].tolist())

    # Store the soil degradation curves in lists
    gg0 = [np.array(back.iloc[0:, i + 1].tolist()) for i in range(len(inp))]

    # Correct degradation curves for negative slope (SIREN)
    shst_norm = [gg0[i] * strn for i in range(len(inp))]
    shst_norm_cor = []
    for i in range(len(shst_norm)):
        temp = shst_norm[i].tolist()
        for j in range(len(temp)):
            if j == 0:
                continue
            elif temp[j - 1] >= temp[j]:
                temp[j] = temp[j - 1] + 1e-7
        temp = np.array(temp)
        shst_norm_cor.append(temp)

    # Define function to calculate damping
    def _get_damping(sh_str_d, strain_lvl_d):
        """
        This function computes hysteretic damping according to the Masing rule.
        
        Parameters
        ----------
        sh_str_d : array
            1D array or list of shear stress {kPa}
        strain_lvl_d : array
            1D array or list of strain levels {%}
            Must be the same length as sh_str_d
            
        Returns
        -------
        damping_d : list
            Hysteretic damping values at the specified strain levels {%}
        """
        from scipy import integrate

        # Obtain the damping curves
        at = sh_str_d * (strain_lvl_d * 0.01) / 2
        int_sh = integrate.cumtrapz(sh_str_d, strain_lvl_d / 100, initial=0)
        al = 8 * (int_sh - at)
        damping_d = al / (4 * np.pi * at)
        for i in range(len(damping_d)):
            if damping_d[i] < 0:
                damping_d[i] = 0
        return damping_d

    # Define final lists of corrected degradation curves
    gg0_cor = [shst_norm_cor[i] / strn for i in range(len(inp))]
    shst_cor = [gg0_cor[i] * inp['G0 (kPa)'][i] * strn * 10
                for i in range(len(inp))]
    dam_cor = [_get_damping(shst_cor[i], strn)
               for i in range(len(inp))]

    return strn, gg0_cor, shst_cor, dam_cor

def _backbone(input_xlsx, out_dir = os.getcwd(), darend_soil='Generic', 
              str_ratio_clay=0.9, str_ratio_sand=0.75, str_haya='Yes', 
              str_vard='Yes', str_match_clay=10.0, str_match_sand=2.0, 
              strain_lvl_dy='', txt=''):
    '''
    This private function derives soil backbone curves following the methods
    from Darendeli (2001) and modifies the peak backbone stresses according to
    Motamed et al. (2016).

    Parameters
    ----------
    input_xlsx : str or path
        Full file path of input spreadsheet (.xlsx)
        See sample on GitHub for formatting
    out_dir : str (optional)
        Directory for results to be written to
        Default = cwd
    darend_soil : str
        Darendeli (2001) soil model
        Options: 'Generic', 'Specific', or 'Generic Sands'
    str_ratio_clay : float (optional)
        Stress to strength ratio for cohesive soils at strain level
        str_match_clay
        Default = 0.9
    str_ratio_sand : float (optional)
        Stress to strength ratio for cohesionless soils at strain level
        str_match_sand
        Default = 0.75
    str_haya : str (optional)
        Option to match stress to strength ratio for cohesionless soils
        after Hayashi (1994)
        Options: 'Yes' or 'No' (default = 'Yes')
    str_vard : str (optional)
        Option to match stress to strength ratio for cohesive soils after
        Vardanega (2012)
        Options: 'Yes' or 'No' (default = 'Yes')
    str_match_clay : float (optional)
        Strain level to match stress to strength ratio for cohesive soils
        Default = 10.0
        Units = %
    str_match_sand : float (optional)
        Strain level to match stress to strength ratio for coheionless
        soils
        Default = 2.0
        Units = %
    strain_lvl_dy : Numpy Array (optional)
        Array of strain values to use for backbone curves
        Values must be float64 with overall shape (10,)
        Units = %
    txt : str (optional)
        Labels for the figures the script produces

    Returns
    -------
    strain_lvl_or : Numpy Array
        Strain levels for backbone data
    gg0_cor : list
        Strength corrected modulus degradation data
    shst_cor : list
        Strength corrected shear stress data
    dam_cor : list
        Masing damping corresponding to gg0_cor and shst_cor
    '''
    
    import init_geo
    from scipy import integrate
    import math
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    
    # Check output directory exists
    if os.path.isdir(out_dir) == False:
        os.makedirs(out_dir)
        
    # Store cwd initialize data object
    params = [input_xlsx, '']
    data = init_geo.update(params)
    depth_wt = data.depth_wt

    # Interpret strain levels
    if strain_lvl_dy == '':
        strain_lvl_dy = np.array([1.e-04, 1.e-03, 3.e-03, 1.e-02, 3.e-02,
                                  1.e-01, 3.e-01, 1.e+00, 3.e+00, 5.e+00])

    # Add darendeli type
    if darend_soil == "Generic_Sands":
        note = "Note: Generalised Darendeli model was used for sands and " +\
               "soil specific for clays"
    elif darend_soil == "Generic":
        note = "Note: Generalised Darendeli model was used for sands and clays"
    elif darend_soil == "Specific":
        note = "Note: Soil specific Darendeli model was used for sands and " +\
               "clays"

    # Define the soil stresses
    mid_dep = data.inp[['Top Depth (m)', 'Bottom Depth (m)']].mean(axis=1)
    tot_str = pd.Series(np.zeros(len(data.inp['Layer #'])), index=data.index)
    uw = pd.Series(np.zeros(len(data.inp['Layer #'])), index=data.index)
    eff_str = pd.Series(np.zeros(len(data.inp['Layer #'])), index=data.index)
    for i in range(data.shape[0]):
        if i == 0:
            tot_str.iloc[i] = mid_dep.iloc[i] * data.inp['gamma (kN/m3)'].iloc[i]
        else:
            tot_str.iloc[i] = tot_str.iloc[i - 1] + \
                data.inp['Thickness (m)'].iloc[i - 1] / 2.0 * \
                data.inp['gamma (kN/m3)'].iloc[i - 1] + \
                data.inp['Thickness (m)'].iloc[i] / 2.0 * \
                data.inp['gamma (kN/m3)'].iloc[i]
    for i in range(data.shape[0]):
        if mid_dep.iloc[i] < depth_wt:
            uw.iloc[i] = 0
        else:
            uw.iloc[i] = (mid_dep.iloc[i] - depth_wt) * 9.81
    eff_str = tot_str - uw
    mean_atm = eff_str * (1 + 2 * data.inp['Ko']) / (3 * 101.3)
    stresses = pd.concat([mid_dep, tot_str, eff_str, mean_atm, uw], axis=1)
    stresses.columns = ['Depth (m)', 'sigma_v (kPa)', 'sigma\'_v (kPa)',
                        'Mean Stress (atm)', 'wat_pr (kPa)']
    # Define dataframe with stress regime
    stress = stresses

    # Define phi coefficients for Darendeli (2001)
    phi1 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi1[i] = 0.0352
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi1[i] = 0.0352
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi1[i] = 0.0416
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi1[i] = 0.0258
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
                    # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi1[i] = 0.0474
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi1[i] = 0.0334
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi1[i] = 0.0416
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi1[i] = 0.0258
    phi2 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi2[i] = 0.00101
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi2[i] = 0.00101
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi2[i] = 0.000689
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi2[i] = 0.00195
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi2[i] = -0.00234
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi2[i] = -0.0000579
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi2[i] = 0.000689
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi2[i] = 0.00195
    phi3 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi3[i] = 0.3246
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi3[i] = 0.3246
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi3[i] = 0.321
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi3[i] = 0.0992
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                  # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi3[i] = 0.25
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi3[i] = 0.249
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi3[i] = 0.321
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi3[i] = 0.0992
    phi4 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi4[i] = 0.3483
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi4[i] = 0.3483
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi4[i] = 0.28
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi4[i] = 0.226
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                    # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi4[i] = 0.234
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi4[i] = 0.482
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi4[i] = 0.28
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi4[i] = 0.226
    phi5 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi5[i] = 0.919
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi5[i] = 0.919
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi5[i] = 1.0
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi5[i] = 0.975
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi5[i] = 0.895
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi5[i] = 0.845
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi5[i] = 1.0
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi5[i] = 0.975
    phi6 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi6[i] = 0.8005
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi6[i] = 0.8005
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi6[i] = 0.712
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi6[i] = 0.958
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                     # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi6[i] = 0.688
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi6[i] = 0.889
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi6[i] = 0.712
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi6[i] = 0.958
    phi7 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi7[i] = 0.0129
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi7[i] = 0.0129
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi7[i] = 0.00303
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi7[i] = 0.00565
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi7[i] = 0.0122
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi7[i] = 0.0202
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi7[i] = 0.00303
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi7[i] = 0.00565
    phi8 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi8[i] = -0.1069
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi8[i] = -0.1069
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi8[i] = -0.1
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi8[i] = -0.1
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                  # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi8[i] = -0.1
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi8[i] = -0.1
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi8[i] = -0.1
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi8[i] = -0.1
    phi9 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi9[i] = -0.2889
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi9[i] = -0.2889
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi9[i] = -0.189
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi9[i] = -0.196
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                     # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi9[i] = -0.127
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi9[i] = -0.372
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi9[i] = -0.189
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi9[i] = -0.196
    phi10 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi10[i] = 0.2919
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi10[i] = 0.2919
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi10[i] = 0.234
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi10[i] = 0.368
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi10[i] = 0.288
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi10[i] = 0.233
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi10[i] = 0.234
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi10[i] = 0.368
    phi11 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi11[i] = 0.6329
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi11[i] = 0.6329
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi11[i] = 0.592
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi11[i] = 0.466
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi11[i] = 0.767
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi11[i] = 0.776
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi11[i] = 0.592
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi11[i] = 0.466
    phi12 = pd.Series(np.zeros(len(range(data.shape[0]))), index=data.index)
    for i in range(data.shape[0]):
        # Darendeli's "Generic" coefficients
        if darend_soil == 'Generic':
            if data.inp['Soil Type'][i] != 'Peat':
                phi12[i] = 0.0057
        # Darendeli's "Generic_Sands" coefficients
        elif darend_soil == 'Generic_Sands':
            if data.inp['Soil Type'][i] != 'Peat':
                if data.inp['Soil Type'][i] == 'Silty Sand'\
                   or data.inp['Soil Type'][i] == 'Sand'\
                   or data.inp['Soil Type'][i] == 'Gravel':
                    phi12[i] = 0.0057
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi12[i] = -0.000767
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi12[i] = 0.0223
        # Darendeli's "Specific" coefficients
        elif darend_soil == 'Specific':
            # Darendeli's "Clean Sands" coefficients
            if data.inp['Soil Type'][i] != 'Peat':
                # Darendeli's "Clean Sands" coefficients
                if data.inp['Soil Type'][i] == 'Sand'\
                        or data.inp['Soil Type'][i] == 'Gravel':
                    phi12[i] = -0.0283
                # Darendeli's 'Sands with high fines content' coefficients
                elif data.inp['Soil Type'][i] == 'Silty Sand':
                    phi12[i] = -0.0294
                # Darendeli's 'Silts' coefficients
                elif data.inp['Soil Type'][i] == 'Silt':
                    phi12[i] = -0.000767
                # Darendeli's 'Clays' coefficients
                elif data.inp['Soil Type'][i] == 'Clay'\
                        or data.inp['Soil Type'][i] == 'Silty Clay':
                    phi12[i] = 0.0223
    # Define parameters to generate G/G0 curves
    beta = pd.Series(np.ones(len(range(data.shape[0]))), index=data.index)
    beta = beta * phi11 + phi12 * np.log(10)
    gamma_r = (phi1 + phi2 * data.inp['PI'] * data.inp['OCR']
               ** phi3) * stress['Mean Stress (atm)']**phi4
    alpha = phi5
    d_min = (phi6 + phi7 * data.inp['PI'] * data.inp['OCR']**phi8) * \
        stress['Mean Stress (atm)']**phi9 * (1 + phi10 * np.log(1))
    daren = pd.concat([stress['Depth (m)'], gamma_r, alpha, d_min], axis=1)
    daren.columns = ['Depth (m)', 'gamma_r', 'alpha', 'Dmin']
    # Obtain the G/G0 curves
    strain_lvl_index = np.array(
        (0.00001,
         0.0001,
         0.0003,
         0.001,
         0.002,
         0.003,
         0.004,
         0.005,
         0.007,
         0.01,
         0.02,
         0.03,
         0.04,
         0.05,
         0.06,
         0.07,
         0.1,
         0.2,
         0.3,
         0.4,
         0.5,
         0.6,
         0.7,
         1,
         2,
         3,
         4,
         5,
         6,
         7,
         10))
    strain_lvl = strain_lvl_index
    # Strains (%)
    g_g0 = [1 / (1 + (strain_lvl / daren['gamma_r'][i])**daren['alpha'][i])
            if data.inp['Soil Type'][i] != 'Peat'
            else np.ones(len(strain_lvl))
            for i in range(daren.shape[0])]
    # Identify peat layers on G/G0 curves
    g_g0_fi = []
    for i in range(data.shape[0]):
        if data.inp['Soil Type'][i] == 'Peat':
            g_g0_fi.append('Peat')
        else:
            g_g0_fi.append(g_g0[i])
    # Calculate the Gmax for each layer (kPa)
    g_max = data.inp['gamma (kN/m3)'] / 9.81 * data.inp['Vs (m/s)']**2
    # Caclulate the shear stress at selected strain levels
    sh_str = [strain_lvl / 100 / (1 + (strain_lvl / daren['gamma_r'][i])) **
              daren['alpha'][i] * g_max[i]
              if data.inp['Soil Type'][i] != 'Peat'
              else np.ones(len(strain_lvl))
              for i in range(daren.shape[0])]
    # Identify peat layers on the shear strength curves
    sh_str_fi = []
    for i in range(data.shape[0]):
        if data.inp['Soil Type'][i] == 'Peat':
            sh_str_fi.append('Peat')
        else:
            sh_str_fi.append(sh_str[i])
    # Obtain the damping curves
    at = [sh_str[i] * (strain_lvl / 100) / 2 for i in range(len(sh_str))]
    int_sh = [integrate.cumtrapz(sh_str[i], strain_lvl / 100, initial=0)
              for i in range(len(sh_str))]
    al = [8 * (int_sh[i] - at[i]) for i in range(len(sh_str))]
    damping = [al[i] / (4 * np.pi * at[i]) + d_min[i] / 100
               if data.inp['Soil Type'][i] != 'Peat'
               else np.ones(len(strain_lvl))
               for i in range(len(sh_str))]
    for damp_curv in damping:
        for i in range(len(damp_curv)):
            if damp_curv[i] < 0:
                damp_curv[i] = 0
    # Identify peat layers on the damping curves
    damping_fi = []
    for i in range(data.shape[0]):
        if data.inp['Soil Type'][i] == 'Peat':
            damping_fi.append('Peat')
        else:
            damping_fi.append(damping[i])
    # Define list of final target curves after Darendeli (2001)
    g_g0_da = g_g0_fi
    sh_str_da = sh_str_fi
    damping_da = damping_fi

    # Hayashi (1994) and Vardanega (2012) target curves
    # Define the strain levels
    strain_lvl = strain_lvl_index
    # Calculate the failure stress (kPa)
    str_fail = np.asarray([data.inp['Su (kPa)'][i] + stress['sigma\'_v (kPa)'][i] *
                           math.tan(math.radians(data.inp['phi (deg)'][i]))
                           for i in range(data.shape[0])])
    # Calculate the Gmax for each layer (kPa)
    g_max = data.inp['gamma (kN/m3)'] / 9.81 * data.inp['Vs (m/s)']**2
    # Define the Hayashi (1994) curves for cohesionless soils
    # Calculate the normalized strain
    gam_gamr = [(strain_lvl / 100) / (str_fail[i] / g_max[i])
                if data.inp['Soil Type'][i] != 'Peat'
                else np.ones(len(strain_lvl))
                for i in range(data.shape[0])]
    # Calculate the shear stress (kPa)
    sh_str_ha = [(np.exp(-0.26 * gam_gamr[i]) * (((2 / 0.3) * gam_gamr[i] + 1)**0.3
                                                 - 1) / (((2 / 0.3) * gam_gamr[i] + 1)**0.3 + 1) + (1 -
                                                                                                    np.exp(-0.26 * gam_gamr[i])) * (((2 * gam_gamr[i] + 1)
                                                                                                                                     - 1) / (2 * gam_gamr[i] + 2))) * str_fail[i]
                 if data.inp['Soil Type'][i] != 'Peat'
                 else np.ones(len(strain_lvl))
                 for i in range(data.shape[0])]
    # Calculate the normalized shear stress
    sh_str_ha_norm = [sh_str_ha[i] / str_fail[i]
                      if data.inp['Soil Type'][i] != 'Peat'
                      else np.ones(len(strain_lvl))
                      for i in range(data.shape[0])]
    # Define the Vardanega (2012) curves for cohesive soils
    # Calculate the strain at which half strength is mobilized
    gamma_m2 = np.asarray(0.004 * (data.inp['OCR'])**0.68)
    # Calculate the coefficient b
    b = np.asarray(0.011 * data.inp['OCR'] + 0.371)
    # Calculate the normalized shear stress
    sh_str_va_norm = [0.5 * ((strain_lvl / 100) / gamma_m2[i])**b[i]
                      for i in range(data.shape[0])]
    # Calculate the normalized shear stress
    sh_str_va = [sh_str_va_norm[i] * str_fail[i]
                 for i in range(data.shape[0])]
    # Add 'nan' values for normalized stress outside the 0.2 to 0.8 range
    for i in range(data.shape[0]):
        for j in range(len(strain_lvl)):
            if 0.2 < sh_str_va_norm[i][j] < 0.8:
                continue
            else:
                sh_str_va_norm[i][j] = np.nan
                sh_str_va[i][j] = np.nan
    # Use curve-fitting to replace 'nan' values in Vardanega (2012)
    sh_str_temp = []
    for i in range(len(sh_str_da)):
        if isinstance(sh_str_da[i], str):
            sh_str_temp.append(np.zeros(len(strain_lvl)))
        else:
            sh_str_temp.append(sh_str_da[i])
    sh_str_in = np.array([sh_str_temp[i][0] for i in range(len(sh_str_temp))])
    sh_str_ta = str_fail
    for i in range(len(sh_str_va)):
        sh_str_va[i][0] = sh_str_in[i]
        sh_str_va[i][-1] = sh_str_ta[i]
    # Get indices with 'nan' values
    nan_ind = [np.argwhere(np.isnan(sh_str_va[i]))
               for i in range(len(sh_str_va))]
    # Get strain levels to apply curve-fitting (currently 'nan' values)
    strain_lvl_cf = [strain_lvl[nan_ind[i]].reshape((-1))
                     for i in range(len(sh_str_va))]
    # Get all indices for strain level vector
    strain_lvl_all = np.arange(len(strain_lvl))
    # Clear 'nan' values in strain data to apply curve fitting
    non_nan_ind = []
    for i in range(len(nan_ind)):
        non_nan_ind.append(np.array([x for x in strain_lvl_all
                                     if x not in nan_ind[i]]))
    strain_lvl_da = [strain_lvl[non_nan_ind[i]].reshape((-1))
                     for i in range(len(sh_str_va))]
    # Clear 'nan' values in shear stress data to apply curve fitting
    sh_str_va_da = [sh_str_va[i][~np.isnan(sh_str_va[i])]
                    for i in range(len(sh_str_va))]
    # Combine the results based on nature of soil
    sh_str = []
    for i in range(data.shape[0]):
        if "Sand" in data.inp['Soil Type'][i] or "Gravel" in data.inp['Soil Type'][i]:
            sh_str.append(sh_str_ha[i])
        elif "Silt" in data.inp['Soil Type'][i] or "Clay" in data.inp['Soil Type'][i]:
            sh_str.append(sh_str_va[i])
        else:
            sh_str.append('Peat')
    # Generate the G/G0 curves
    g_g0 = []
    for i in range(data.shape[0]):
        if isinstance(sh_str[i], str):
            g_g0.append('Peat')
        else:
            g_g0.append(sh_str[i] / (strain_lvl / 100) / g_max[i])
    # Define list of final target curves after Hayashi (1994)
    # and Vardanega (2012)
    g_g0_hv = g_g0
    sh_str_hv = sh_str

    # Derive curves for organic soil according to Kishida (1996) and Kishida
    # et al. (2009)

    # Define input parameters required to derive degradation curve
    oc = np.asarray(data.inp['OC (%)'].astype(float))
    vs = np.asarray(data.inp['Vs (m/s)'].astype(float))
    g_max = np.asarray(data.inp['gamma (kN/m3)'] / 9.81 * data.inp['Vs (m/s)']**2)
    OCR = np.asarray(data.inp['OCR'])
    LCR = OCR
    stre = np.asarray(stress['sigma\'_v (kPa)'])
    su = np.asarray(data.inp['Su (kPa)'])
    # Define the strain levels
    strain_lvl_or = strain_lvl_dy  # Strains (%)
    # Assess dergadation curve based on Kishida et al. (2006) & (2009)
    g_g0 = []
    for i in range(data.shape[0]):
        if data.inp['Soil Type'][i] != 'Peat':
            g_g0.append('No Peat')
        else:
            lna = 5.2 + 0.48 * (2 / (1 + math.exp(oc[i] / 23))) + 0.74 *\
                ((((6 / (1 + math.exp(oc[i] / 23))) - 1.5)
                  / math.log(1 + 3 * math.exp(1 +
                                              (6 / (1 + math.exp(oc[i] / 23)))))) - 1)
            m = 0.8 - 0.4 * (2 / (1 + math.exp(oc[i] / 23)))
            n = 1 - 0.37 * (2 / (1 + math.exp(oc[i] / 23)))
            X2 = math.log(stre[i])
            X3 = 2 / (1 + math.exp(oc[i] / 23))
            X4 = math.log(LCR[i])
            b4, b0, b1, b3, b6, = 0.8 - 0.4 * X3, 5.110, -0.729, -0.693, 0.0
            b9, b10 = -1.41, -0.95
            ref_str = math.exp(b9 + b10 * (X3 - 0.5))
            X1 = np.log(strain_lvl_or + ref_str)
            b2 = 1 - 0.185 * (1 + ((np.log(ref_str) + 2.5) /
                                   np.log((1 + strain_lvl_or) / ref_str)))
            b5 = 0.185 / np.log((1 + strain_lvl_or) / ref_str)
            b7 = -0.37 * (1 + ((np.log(ref_str) + 2.5)) /
                          np.log((1 + strain_lvl_or) / ref_str))
            b8 = 0.37 / np.log((1 + strain_lvl_or) / ref_str)
            lng = b0 + b1 * X1 + b2 * X2 + b3 * X3 + b4 * X4 + b5 * (X2 - 4) * (X3 - 0.5) +\
                b6 * (X1 + 2.5) * (X3 - 0.5)\
                + b7 * (X2 - 4) * (X3 - 0.5) + b8 * \
                (X1 + 2.5) * (X2 - 4) * (X3 - 0.5)
            thg = np.exp(lng)
            g_g0_or = thg / np.max(thg)
            g_g0.append(g_g0_or)
    # Define shear stress strain curves for the organic layers
    sh_str = []
    for i in range(data.shape[0]):
        if data.inp['Soil Type'][i] != 'Peat':
            sh_str.append('No Peat')
        else:
            sh_str.append(g_max[i] * strain_lvl_or)
    # Define list of final target curves after Kishida (1996)
    g_g0_or = g_g0
    sh_str_or = sh_str

    # Fit hyperbola to target curves
    # Define list with all target curves
    g_g0_da_all = g_g0_da
    sh_str_da_all = sh_str_da
    damping_da_all = damping_da
    g_g0_hv_all = g_g0_hv
    sh_str_hv_all = sh_str_hv

    # Create an empty list to store results
    gg0 = []
    shst = []
    damp = []

    # Define function to minimize residuals between Groholski (2016) and the
    # target curves for small and large strains
    def objective(x, curve_id, strain_lvl, stress, data, str_bn_fit,
                  str_bn_fit_2, g_g0_da, g_g0_hv, soil_type, ind_val=0):

        # Initialize Groholski (2016) fitting parameters
        theta_1 = x[0]
        theta_2 = x[1]
        theta_3 = x[2]
        # Calculate normalized stress-strain curve by Groholski (2016)
        tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] * math.tan(math.radians(
            data.inp['phi (deg)']. iloc[curve_id])) + data.inp['Su (kPa)'].iloc[curve_id]
        G_max = data.inp['G0 (kPa)'].iloc[curve_id]
        gamma_r = tau_max / G_max
        theta_tau = theta_1 + theta_2 * \
            (strain_lvl / 100 / gamma_r) / \
            (theta_3 + (strain_lvl / 100 / gamma_r))
        theta_tau = np.asarray(theta_tau)
        theta_tau = np.minimum(theta_tau, 1)
        sh_str_ha_no = 2 * (strain_lvl / 100 / gamma_r) / (1 +
                                                           (strain_lvl / 100 / gamma_r) +
                                                           ((1 + (strain_lvl / 100 / gamma_r))**2 -
                                                            4 * theta_tau * (strain_lvl / 100 / gamma_r))**0.5)
        sh_str_ha = sh_str_ha_no * tau_max
        g_g0_ha = sh_str_ha / (strain_lvl / 100) / G_max
        # Calculate residulas up to target strain
        ind_1 = int(np.where(strain_lvl == str_bn_fit)[0])
        # Obtain residulas for Darendeli (2001) region
        res_1 = np.sum(np.square(g_g0_da[0:ind_1] - g_g0_ha[0:ind_1]))
        if soil_type == 'cohesive':
            # Calculate residulas for Vardanega points
            g_g0_ha_red = g_g0_ha[ind_val]
            res_2 = np.sum(np.square(g_g0_hv - g_g0_ha_red))
        else:
            # Obtain residulas for Hayashi region
            ind_2 = int(np.where(strain_lvl == str_bn_fit)[0])
            res_2 = np.sum(np.square(g_g0_hv[ind_2:] - g_g0_ha[ind_2:]))
        # Define weighting coefficients on target curves
        a1 = 1.0
        a2 = 1.0
        # Define residulas to be minimized
        res = a1 * res_1 + a2 * res_2

        return res

    def cons1(x):
        '''Constant defined as: theta_1 + theta_2 < 1.0'''
        theta_1 = x[0]
        theta_2 = x[1]
        return 1 - theta_1 - theta_2

    def cons2(x, curve_id, strain_lvl, stress, data, str_bn_fit,
              str_bn_fit_2, g_g0_da):
        '''Constant defined to minimize distance between
        Darendeli (2001) and Groholski (2016) curves. Rule is the
        average vertical spacing between the two curves up to
        maximum defined strain level to be equal to mean_spac'''
        # Define this value if you want to change constrain
        mean_spac = 0.001
        # Initialize Groholski (2016) fitting parameters
        theta_1 = x[0]
        theta_2 = x[1]
        theta_3 = x[2]
        # Calculate normalized stress-strain curve by Groholski (2016)
        tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] * math.tan(math.radians(
            data.inp['phi (deg)']. iloc[curve_id])) + data.inp['Su (kPa)'].iloc[curve_id]
        G_max = data.inp['G0 (kPa)'].iloc[curve_id]
        gamma_r = tau_max / G_max
        theta_tau = theta_1 + theta_2 * \
            (strain_lvl / 100 / gamma_r) / \
            (theta_3 + (strain_lvl / 100 / gamma_r))
        theta_tau = np.asarray(theta_tau)
        theta_tau = np.minimum(theta_tau, 1)
        sh_str_ha_no = 2 * (strain_lvl / 100 / gamma_r) / (1 +
                                                           (strain_lvl / 100 / gamma_r) +
                                                           ((1 + (strain_lvl / 100 / gamma_r))**2 -
                                                            4 * theta_tau * (strain_lvl / 100 / gamma_r))**0.5)
        sh_str_ha = sh_str_ha_no * tau_max
        g_g0_ha = sh_str_ha / (strain_lvl / 100) / G_max
        # Calculate residuals up to target strain
        ind = int(np.where(strain_lvl == str_bn_fit)[0])
        # Defined constrain rule to be met during optimization
        diff = mean_spac / len(g_g0_da[:ind]) -\
            np.sum(np.abs(g_g0_ha[:ind] - g_g0_da[:ind]))
        return diff

    def cons3(x, curve_id, strain_lvl, stress, data, str_bn_fit,
              str_bn_fit_2, g_g0_hv, ind_diff, str_ratio):
        '''Constant to limit the maximum vertical distance on the
        stress-strain curve, between the shear strength and the stress
        at a given strain level. Alternatively, you can set this
        constrain to minimize vertical distance between Hayash (1994)
        and Vardanega (2012) to Darendeli (2001)'''
        # Define this value if you want to change constrain
        #mean_spac = 0.01
        # Initialize Groholski (2016) fitting parameters
        theta_1 = x[0]
        theta_2 = x[1]
        theta_3 = x[2]
        # Calculate normalized stress-strain curve by Groholski (2016)
        tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] * math.tan(math.radians(
            data.inp['phi (deg)']. iloc[curve_id])) + data.inp['Su (kPa)'].iloc[curve_id]
        G_max = data.inp['G0 (kPa)'].iloc[curve_id]
        gamma_r = tau_max / G_max
        theta_tau = theta_1 + theta_2 * \
            (strain_lvl / 100 / gamma_r) / \
            (theta_3 + (strain_lvl / 100 / gamma_r))
        theta_tau = np.asarray(theta_tau)
        theta_tau = np.minimum(theta_tau, 1)
        sh_str_ha_no = 2 * (strain_lvl / 100 / gamma_r) / (1 +
                                                           (strain_lvl / 100 / gamma_r) +
                                                           ((1 + (strain_lvl / 100 / gamma_r))**2 -
                                                            4 * theta_tau * (strain_lvl / 100 / gamma_r))**0.5)
        diff = sh_str_ha_no[ind_diff] - str_ratio
        return diff

    # Loop through all curves to fit hyperbolic Groholski (2016) model
    for k in range(len(data)):

        # Import curve to be examined
        curve_id = int(k)
        curve_id_print = int(k) + 1
        print('Soil layer: ' + str(curve_id_print))

        if 'Peat' in data.inp['Soil Type'].iloc[curve_id]:
            gg0.append("Peat")
            shst.append("Peat")
            damp.append("Peat")

        elif 'Clay' in data.inp['Soil Type'].iloc[curve_id]\
                or 'Silt' == data.inp['Soil Type'].iloc[curve_id]:
            # Define target curves for individual layers
            g_g0_da = g_g0_da_all[curve_id]
            sh_str_da = sh_str_da_all[curve_id]
            damping_da = damping_da_all[curve_id]
            g_g0_hv = g_g0_hv_all[curve_id]
            sh_str_hv = sh_str_hv_all[curve_id]
            # Define the strain level that peak strength is reached
            if str_vard == 'Yes':
                tau_max = data.inp['Su (kPa)'].iloc[curve_id]
                diff = np.abs(sh_str_hv - tau_max) / tau_max * 100
                for i in range(len(diff)):
                    if (math.isnan(diff[i])):
                        diff[i] = 1e20
                ind_diff = np.where(diff < 5)
                strain_cut = strain_lvl[ind_diff][0]
                ind_diff = np.where(strain_lvl == strain_cut)
                ind_diff = int(ind_diff[0])
            elif str_vard == 'No':
                ind_clay = np.argmin(np.abs(strain_lvl - str_match_clay))
                ind_diff = ind_clay
            # Get indices with nan values on strains
            ind_nan = np.argwhere(np.isnan(g_g0_hv))
            # Get strains with nan values
            strn_nan = strain_lvl[ind_nan]
            strn_nan = [item for sublist in strn_nan for item in sublist]
            ind_nan = [item for sublist in ind_nan for item in sublist]
            ind_val = sorted(set(range(len(strain_lvl))) - set(ind_nan))
            # Get shress strain curves without nan values
            set_diff = sorted(set(strain_lvl) - set(strn_nan))
            strn_val = np.array(set_diff)
            strn_nan = np.array(strn_nan)
            # Remove nan values from Vardanega
            g_g0_hv = g_g0_hv[np.logical_not(np.isnan(g_g0_hv))]
            sh_str_hv = sh_str_hv[np.logical_not(np.isnan(sh_str_hv))]
            # Initialize Groholski (2016) fitting parameters
            theta_1 = 0.5
            theta_2 = 0.5
            theta_3 = 0.5
            # Optimize curves to fit target curves at small strain
            str_bn_fit = 0.1  # Define boundary of max strain for target curve
            str_bn_fit_2 = 0.1  # Define boundary of max strain for target curve
            x = np.array([theta_1, theta_2, theta_3])  # Define variable array

            # Run optimization
            # Initialize Groholski (2016) fitting parameters
            theta_1 = 0.5
            theta_2 = 0.5
            theta_3 = 0.5
            x0 = np.array([theta_1, theta_2, theta_3])
            # Set the constrains in a dictionary format
            con1 = {'type': 'ineq', 'fun': cons1}
            con2 = {'type': 'ineq', 'fun': cons2,
                    'args': (curve_id, strain_lvl, stress, data, str_bn_fit,
                             str_bn_fit_2, g_g0_da,)}
            con3 = {'type': 'ineq', 'fun': cons3,
                    'args': (curve_id, strain_lvl, stress, data, str_bn_fit,
                             str_bn_fit_2, g_g0_hv, ind_diff, str_ratio_clay)}
            cons = ([con1, con2, con3])
            # Obtain the solution
            solution = minimize(
                objective,
                x0,
                args=(
                    curve_id,
                    strain_lvl,
                    stress,
                    data,
                    str_bn_fit,
                    str_bn_fit_2,
                    g_g0_da,
                    g_g0_hv,
                    'cohesive',
                    ind_val),
                method='SLSQP',
                constraints=cons,
                options={
                    'disp': True,
                    'maxiter': 2000})
            # Obtain the fitting parameters
            x = solution.x
            # Get final curves after Groholski (2016)
            theta_1 = x[0]
            theta_2 = x[1]
            theta_3 = x[2]
            tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] * math.tan(math.radians(
                data.inp['phi (deg)']. iloc[curve_id])) + data.inp['Su (kPa)'].iloc[curve_id]
            G_max = data.inp['G0 (kPa)'].iloc[curve_id]
            gamma_r = tau_max / G_max
            theta_tau = theta_1 + theta_2 * \
                (strain_lvl / 100 / gamma_r) / \
                (theta_3 + (strain_lvl / 100 / gamma_r))
            theta_tau = np.asarray(theta_tau)
            theta_tau = np.minimum(theta_tau, 1)
            sh_str_ha_no = 2 * (strain_lvl / 100 / gamma_r) / (1 +
                                                               (strain_lvl / 100 / gamma_r) +
                                                               ((1 + (strain_lvl / 100 / gamma_r))**2 -
                                                                4 * theta_tau * (strain_lvl / 100 / gamma_r))**0.5)
            sh_str_ha = sh_str_ha_no * tau_max
            g_g0_ha = sh_str_ha / (strain_lvl / 100) / G_max
            # Obtain the Masing damping curve for Groholski (2016) curve
            at = sh_str_ha * (strain_lvl / 100) / 2
            int_sh = integrate.cumtrapz(sh_str_ha, strain_lvl / 100, initial=0)
            al = 8 * (int_sh - at)
            damping = al / (4 * np.pi * at)
            # Remove negative values from damping curve
            for i in range(len(damping)):
                if damping[i] < 0:
                    damping[i] = 0
            damping_ha = damping

            # Check for non-monotonically decreasing degradation curve
            while not np.all(np.diff(g_g0_ha) < 0) or\
                    not np.all(np.diff(sh_str_ha) > 0):
                str_ratio_clay = str_ratio_clay - 0.05
                # Initialize Groholski (2016) fitting parameters
                theta_1 = 0.5
                theta_2 = 0.5
                theta_3 = 0.5
                x0 = np.array([theta_1, theta_2, theta_3])
                # Set the constrains in a dictionary format
                con1 = {'type': 'ineq', 'fun': cons1}
                con2 = {
                    'type': 'ineq',
                    'fun': cons2,
                    'args': (
                        curve_id,
                        strain_lvl,
                        stress,
                        data,
                        str_bn_fit,
                        str_bn_fit_2,
                        g_g0_da,
                    )}
                con3 = {
                    'type': 'ineq',
                    'fun': cons3,
                    'args': (
                        curve_id,
                        strain_lvl,
                        stress,
                        data,
                        str_bn_fit,
                        str_bn_fit_2,
                        g_g0_hv,
                        ind_diff,
                        str_ratio_clay)}
                cons = ([con1, con2, con3])
                # Obtain the solution
                solution = minimize(
                    objective,
                    x0,
                    args=(
                        curve_id,
                        strain_lvl,
                        stress,
                        data,
                        str_bn_fit,
                        str_bn_fit_2,
                        g_g0_da,
                        g_g0_hv,
                        'cohesive',
                        ind_val),
                    method='SLSQP',
                    constraints=cons,
                    options={
                        'disp': True,
                        'maxiter': 2000})
                # Obtain the fitting parameters
                x = solution.x
                # Get final curves after Groholski (2016)
                theta_1 = x[0]
                theta_2 = x[1]
                theta_3 = x[2]
                tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] * math.tan(math.radians(
                    data.inp['phi (deg)']. iloc[curve_id])) + data.inp['Su (kPa)'].iloc[curve_id]
                G_max = data.inp['G0 (kPa)'].iloc[curve_id]
                gamma_r = tau_max / G_max
                theta_tau = theta_1 + theta_2 * \
                    (strain_lvl / 100 / gamma_r) / \
                    (theta_3 + (strain_lvl / 100 / gamma_r))
                theta_tau = np.asarray(theta_tau)
                theta_tau = np.minimum(theta_tau, 1)
                sh_str_ha_no = 2 * (strain_lvl / 100 / gamma_r) / (1 +
                                                                   (strain_lvl / 100 / gamma_r) +
                                                                   ((1 + (strain_lvl / 100 / gamma_r))**2 -
                                                                    4 * theta_tau * (strain_lvl / 100 / gamma_r))**0.5)
                sh_str_ha = sh_str_ha_no * tau_max
                g_g0_ha = sh_str_ha / (strain_lvl / 100) / G_max

                # Obtain the Masing damping curve for Groholski (2016) curve
                at = sh_str_ha * (strain_lvl / 100) / 2
                int_sh = integrate.cumtrapz(
                    sh_str_ha, strain_lvl / 100, initial=0)
                al = 8 * (int_sh - at)
                damping = al / (4 * np.pi * at)
                # Remove negative values from damping curve
                for i in range(len(damping)):
                    if damping[i] < 0:
                        damping[i] = 0
                damping_ha = damping

            # Plot results for cohesive soils
            # Create a figure instance
            plt.figure

            # Define axes instances
            ax1 = plt.subplot2grid((82, 116), (7, 2), colspan=51, rowspan=29)
            ax2 = plt.subplot2grid((82, 116), (7, 65), colspan=51, rowspan=29)
            ax3 = plt.subplot2grid((82, 116), (46, 2), colspan=51, rowspan=29)
            ax4 = plt.subplot2grid((82, 116), (46, 65), colspan=51, rowspan=29)

            # Plot degradation curves
            ax1.semilogx(strain_lvl, g_g0_da, color='r',
                         marker='o', label='Darendeli')
            ax1.semilogx(strain_lvl, g_g0_ha, color='b',
                         marker='o', label='Groholski')
            ax1.set_ylim([0.0, 1.0])
            ax1.set_xlim([1e-5, 10.0])
            ax1.set_xlabel('Shear strain (%)')
            ax1.set_ylabel('G/G0')
            ax1.legend(loc='lower left')
            ax1.set_title('Shear modulus degradation curve')

            # Plot stress-strain curves
            ax2.plot(strain_lvl, sh_str_da, color='r',
                     marker='o', label='Darendeli')
            ax2.plot(strain_lvl, sh_str_ha, color='b',
                     marker='o', label='Groholski')
            ax2.plot(strn_val, sh_str_hv, color='g',
                     marker='o', label='Vardanega')
            ax2.set_xlabel('Shear strain (%)')
            ax2.set_ylabel('Shear stress (kPa)')
            ax2.set_xlim([0.0, 10.0])
            ax2.set_ylim(bottom=0)
            ax2.legend(loc='lower right')
            ax2.set_title('Backbone curve')

            # Plot damping curves
            ax3.semilogx(strain_lvl, damping_da * 100, color='r',
                         marker='o', label='Darendeli')
            ax3.semilogx(strain_lvl, damping_ha * 100, color='b',
                         marker='o', label='Groholski')
            ax3.set_xlim([1e-5, 10.0])
            ax3.set_ylim([0.0, 100.0])
            ax3.set_xlabel('Shear strain (%)')
            ax3.set_ylabel('Damping (%)')
            ax3.legend(loc='upper left')
            ax3.set_title('Damping curve')

            # Plot theta_tau curves
            ax4.plot(strain_lvl, theta_tau, color='b',
                     marker='o', label='Groholski')
            ax4.set_xlabel('Shear strain (%)')
            ax4.set_ylabel('Theta tau (-)')
            ax4.set_xlim([0.0, 10.0])
            ax4.legend(loc='upper right')
            ax4.set_title('Theta-tau curve')

            # Add plot title
            flarge = round(1.50 * plt.rcParams['font.size'])
            fsmall = round(0.85 * plt.rcParams['font.size'])
            title = "Soil Layer " + str(k + 1) + " - " + data.inp['Soil Type'][k]
            plt.suptitle(title, fontsize=flarge, fontweight='bold',
                         horizontalalignment='left',
                         verticalalignment='top', x=0.02, y=0.97)

            plt.figtext(0.98, 0.97, txt, fontsize=fsmall,
                        horizontalalignment='right',
                        verticalalignment='top')

            plt.figtext(0.02, 0.03, note, fontstyle='oblique',
                        horizontalalignment='left',
                        verticalalignment='bottom')

            # Save figures
            name = os.path.join(out_dir, title)
            plt.savefig(name)
            plt.close()

            # Store final curves
            gg0.append(g_g0_ha)
            shst.append(sh_str_ha)
            damp.append(damping_ha)

        # Repeat for cohesionless soils
        elif 'Sand' in data.inp['Soil Type'].iloc[curve_id]\
                or 'Gravel' == data.inp['Soil Type'].iloc[curve_id]:

            # Define target curves for individual layers
            g_g0_da = g_g0_da_all[curve_id]
            sh_str_da = sh_str_da_all[curve_id]
            damping_da = damping_da_all[curve_id]
            g_g0_hv = g_g0_hv_all[curve_id]
            sh_str_hv = sh_str_hv_all[curve_id]

            # Define the strain level that peak strength is reached
            if str_haya == 'Yes':
                tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] *\
                    math.tan(math.radians(data.inp['phi (deg)'].iloc[curve_id]))\
                    + data.inp['Su (kPa)'].iloc[curve_id]
                diff = np.abs(sh_str_hv - tau_max) / tau_max * 100
                for i in range(len(diff)):
                    if (math.isnan(diff[i])):
                        diff[i] = 1e20
                ind_diff = np.where(diff < 5)
                strain_cut = strain_lvl[ind_diff][0]
                ind_diff = np.where(strain_lvl == strain_cut)
                ind_diff = int(ind_diff[0])
            elif str_haya == 'No':
                ind_sand = np.argmin(np.abs(strain_lvl - str_match_sand))
                ind_diff = ind_sand

            # Initialize Groholski (2016) fitting parameters
            theta_1 = 0.5
            theta_2 = 0.5
            theta_3 = 0.5
            # Optimize curves to fit target curves at small strain
            str_bn_fit = 0.1  # Define boundary of max strain for target curve
            str_bn_fit_2 = 0.1  # Define boundary of max strain for target curve
            x = np.array([theta_1, theta_2, theta_3])  # Define variable array

            # Run optimization
            # Initialize Groholski (2016) fitting parameters
            theta_1 = 0.5
            theta_2 = 0.5
            theta_3 = 0.5
            x0 = np.array([theta_1, theta_2, theta_3])
            # Set the constrains in a dictionary format
            con1 = {'type': 'ineq', 'fun': cons1}
            con2 = {'type': 'ineq', 'fun': cons2,
                    'args': (curve_id, strain_lvl, stress, data, str_bn_fit,
                             str_bn_fit_2, g_g0_da,)}
            con3 = {'type': 'ineq', 'fun': cons3,
                    'args': (curve_id, strain_lvl, stress, data, str_bn_fit,
                             str_bn_fit_2, g_g0_hv, ind_diff, str_ratio_sand)}
            cons = ([con1, con2, con3])
            # Obtain the solution
            solution = minimize(
                objective,
                x0,
                args=(
                    curve_id,
                    strain_lvl,
                    stress,
                    data,
                    str_bn_fit,
                    str_bn_fit_2,
                    g_g0_da,
                    g_g0_hv,
                    'cohesionless'),
                method='SLSQP',
                constraints=cons,
                options={
                    'disp': True,
                    'maxiter': 2000})
            # Obtain the fitting parameters
            x = solution.x
            # Get final curves after Groholski (2016)
            theta_1 = x[0]
            theta_2 = x[1]
            theta_3 = x[2]
            tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] * math.tan(math.radians(
                data.inp['phi (deg)']. iloc[curve_id])) + data.inp['Su (kPa)'].iloc[curve_id]
            G_max = data.inp['G0 (kPa)'].iloc[curve_id]
            gamma_r = tau_max / G_max
            theta_tau = theta_1 + theta_2 * \
                (strain_lvl / 100 / gamma_r) / \
                (theta_3 + (strain_lvl / 100 / gamma_r))
            theta_tau = np.asarray(theta_tau)
            theta_tau = np.minimum(theta_tau, 1)
            sh_str_ha_no = 2 * (strain_lvl / 100 / gamma_r) / (1 +
                                                               (strain_lvl / 100 / gamma_r) +
                                                               ((1 + (strain_lvl / 100 / gamma_r))**2 -
                                                                4 * theta_tau * (strain_lvl / 100 / gamma_r))**0.5)
            sh_str_ha = sh_str_ha_no * tau_max
            g_g0_ha = sh_str_ha / (strain_lvl / 100) / G_max
            # Obtain the Masing damping curve for Groholski (2016) curve
            at = sh_str_ha * (strain_lvl / 100) / 2
            int_sh = integrate.cumtrapz(sh_str_ha, strain_lvl / 100, initial=0)
            al = 8 * (int_sh - at)
            damping = al / (4 * np.pi * at)
            # Remove negative values from damping curve
            for i in range(len(damping)):
                if damping[i] < 0:
                    damping[i] = 0
            damping_ha = damping

            # Check for non-monotonically decreasing degradation curve
            while not np.all(np.diff(g_g0_ha) < 0) or\
                    not np.all(np.diff(sh_str_ha) > 0):
                str_ratio_sand = str_ratio_sand - 0.05
                # Initialize Groholski (2016) fitting parameters
                theta_1 = 0.5
                theta_2 = 0.5
                theta_3 = 0.5
                x0 = np.array([theta_1, theta_2, theta_3])
                # Set the constrains in a dictionary format
                con1 = {'type': 'ineq', 'fun': cons1}
                con2 = {
                    'type': 'ineq',
                    'fun': cons2,
                    'args': (
                        curve_id,
                        strain_lvl,
                        stress,
                        data,
                        str_bn_fit,
                        str_bn_fit_2,
                        g_g0_da,
                    )}
                con3 = {
                    'type': 'ineq',
                    'fun': cons3,
                    'args': (
                        curve_id,
                        strain_lvl,
                        stress,
                        data,
                        str_bn_fit,
                        str_bn_fit_2,
                        g_g0_hv,
                        ind_diff,
                        str_ratio_sand)}
                cons = ([con1, con2, con3])
                # Obtain the solution
                solution = minimize(
                    objective,
                    x0,
                    args=(
                        curve_id,
                        strain_lvl,
                        stress,
                        data,
                        str_bn_fit,
                        str_bn_fit_2,
                        g_g0_da,
                        g_g0_hv,
                        'cohesionless'),
                    method='SLSQP',
                    constraints=cons,
                    options={
                        'disp': True,
                        'maxiter': 2000})
                # Obtain the fitting parameters
                x = solution.x
                # Get final curves after Groholski (2016)
                theta_1 = x[0]
                theta_2 = x[1]
                theta_3 = x[2]
                tau_max = stress['sigma\'_v (kPa)'].iloc[curve_id] * math.tan(math.radians(
                    data.inp['phi (deg)']. iloc[curve_id])) + data.inp['Su (kPa)'].iloc[curve_id]
                G_max = data.inp['G0 (kPa)'].iloc[curve_id]
                gamma_r = tau_max / G_max
                theta_tau = theta_1 + theta_2 * \
                    (strain_lvl / 100 / gamma_r) / \
                    (theta_3 + (strain_lvl / 100 / gamma_r))
                theta_tau = np.asarray(theta_tau)
                theta_tau = np.minimum(theta_tau, 1)
                sh_str_ha_no = 2 * (strain_lvl / 100 / gamma_r) / (1 +
                                                                   (strain_lvl / 100 / gamma_r) +
                                                                   ((1 + (strain_lvl / 100 / gamma_r))**2 -
                                                                    4 * theta_tau * (strain_lvl / 100 / gamma_r))**0.5)
                sh_str_ha = sh_str_ha_no * tau_max
                g_g0_ha = sh_str_ha / (strain_lvl / 100) / G_max

                # Obtain the Masing damping curve for Groholski (2016) curve
                at = sh_str_ha * (strain_lvl / 100) / 2
                int_sh = integrate.cumtrapz(
                    sh_str_ha, strain_lvl / 100, initial=0)
                al = 8 * (int_sh - at)
                damping = al / (4 * np.pi * at)
                # Remove negative values from damping curve
                for i in range(len(damping)):
                    if damping[i] < 0:
                        damping[i] = 0
                damping_ha = damping

            # Plot defined and target curves
            # Create a figure instance
            plt.figure

            # Define axes instances
            ax1 = plt.subplot2grid((82, 116), (7, 2), colspan=51, rowspan=29)
            ax2 = plt.subplot2grid((82, 116), (7, 65), colspan=51, rowspan=29)
            ax3 = plt.subplot2grid((82, 116), (46, 2), colspan=51, rowspan=29)
            ax4 = plt.subplot2grid((82, 116), (46, 65), colspan=51, rowspan=29)

            # Plot degradation curves
            ax1.semilogx(strain_lvl, g_g0_da, color='r',
                         marker='o', label='Darendeli')
            ax1.semilogx(strain_lvl, g_g0_ha, color='b',
                         marker='o', label='Groholski')
            ax1.set_xlim([1e-5, 10.0])
            ax1.set_ylim([0.0, 1.0])
            ax1.set_xlabel('Shear strain (%)')
            ax1.set_ylabel('G/G0')
            ax1.legend(loc='lower left')
            ax1.set_title('Shear modulus degradation curve')

            # Plot stress-strain curves
            ax2.plot(strain_lvl, sh_str_da, color='r',
                     marker='o', label='Darendeli')
            ax2.plot(strain_lvl, sh_str_ha, color='b',
                     marker='o', label='Groholski')
            ax2.plot(strain_lvl, sh_str_hv, color='g',
                     marker='o', label='Hayashi')
            ax2.set_xlabel('Shear strain (%)')
            ax2.set_ylabel('Shear stress (kPa)')
            ax2.set_xlim([0.0, 10.0])
            ax2.set_ylim(bottom=0)
            ax2.legend(loc='lower right')
            ax2.set_title('Backbone curve')

            # Plot damping curves
            ax3.semilogx(strain_lvl, damping_da * 100, color='r',
                         marker='o', label='Darendeli')
            ax3.semilogx(strain_lvl, damping_ha * 100, color='b',
                         marker='o', label='Groholski')
            ax3.set_xlim([1e-5, 10.0])
            ax3.set_ylim([0.0, 100.0])
            ax3.set_xlabel('Shear strain (%)')
            ax3.set_ylabel('Damping (%)')
            ax3.legend(loc='upper left')
            ax3.set_title('Damping curve')

            # Plot theta_tau curves
            ax4.plot(strain_lvl, theta_tau, color='b',
                     marker='o', label='Groholski')
            ax4.set_xlabel('Shear strain (%)')
            ax4.set_ylabel('Theta tau (-)')
            ax4.set_ylim([0.0, 1.0])
            ax4.legend(loc='upper right')
            ax4.set_title('Theta-tau curve')

            # Add plot title
            flarge = round(1.50 * plt.rcParams['font.size'])
            fsmall = round(0.85 * plt.rcParams['font.size'])

            title = "Soil Layer " + str(k + 1) + " - " + data.inp['Soil Type'][k]
            plt.suptitle(title, fontsize=flarge, fontweight='bold',
                         horizontalalignment='left',
                         verticalalignment='top', x=0.02, y=0.97)

            plt.figtext(0.98, 0.97, txt, fontsize=fsmall,
                        horizontalalignment='right',
                        verticalalignment='top')

            plt.figtext(0.02, 0.03, note, fontstyle='oblique',
                        horizontalalignment='left',
                        verticalalignment='bottom')

            # Save figures
            name = os.path.join(out_dir, title)
            plt.savefig(name)
            plt.close()

            # Store final curves
            gg0.append(g_g0_ha)
            shst.append(sh_str_ha)
            damp.append(damping_ha)

    # Define derived curves for all layers
    g_g0_ha = gg0
    sh_str_ha = shst
    damping_ha = damp
    g_g0_or = g_g0_or
    sh_str_or = sh_str_or

    # Define function to calculate damping
    def get_damping(sh_str_d, strain_lvl_d):
        """
        This function computes Masing hysteretic damping given
        the stress-strain backbone curve.
        """
        # Import libraries
        from scipy import integrate
        # Obtain the damping curves
        at = sh_str_d * (strain_lvl_d * 0.01) / 2
        int_sh = integrate.cumtrapz(sh_str_d, strain_lvl_d / 100, initial=0)
        al = 8 * (int_sh - at)
        damping_d = al / (4 * np.pi * at)
        for i in range(len(damping_d)):
            if damping_d[i] < 0:
                damping_d[i] = 0
        return damping_d

    # Define number of layers
    num_layers = len(data)

    # Adjust strain level to non-organic layers
    index = np.argsort(strain_lvl)
    sorted_x = strain_lvl[index]
    sorted_index = np.searchsorted(sorted_x, strain_lvl_or)
    yindex = np.take(index, sorted_index, mode="clip")
    mask = strain_lvl[yindex] != strain_lvl_or
    index = np.ma.array(yindex, mask=mask)
    index = [index[i] for i in range(len(index))]

    # Define empty lists to store final curves
    gg0 = []
    shst = []
    damp = []

    # Replace data for peat layers
    for i in range(num_layers):
        if data.inp['Soil Type'][i] == 'Peat':
            gg0.append(g_g0_or[i])
            shst.append(sh_str_or[i])
        else:
            gg0.append(g_g0_ha[i][index])
            shst.append(sh_str_ha[i][index])
            damp.append(damping_ha[i][index])

    # Correct degradation curves for negative slope (SIREN)
    strn = strain_lvl_or
    shst_norm = [gg0[i] * strn for i in range(num_layers)]
    shst_norm_cor = []
    for i in range(len(shst_norm)):
        temp = shst_norm[i].tolist()
        for j in range(len(temp)):
            if j == 0:
                continue
            elif temp[j - 1] >= temp[j]:
                temp[j] = temp[j - 1] + 1e-7
        temp = np.array(temp)
        shst_norm_cor.append(temp)

    # Define final lists of corrected degradation curves
    gg0_cor = [shst_norm_cor[i] / strn for i in range(num_layers)]
    shst_cor = [gg0_cor[i] * data.inp['G0 (kPa)'][i] * strn * 10
                for i in range(num_layers)]
    dam_cor = [get_damping(shst_cor[i], strn)
               for i in range(num_layers)]

    # Plot results for peat layers
    for k in range(data.shape[0]):
        if data.inp['Soil Type'][k] == 'Peat':
            # Create a figure instance
            plt.figure

            # Define axes instances
            ax1 = plt.subplot2grid((82, 116), (7, 2), colspan=51, rowspan=29)
            ax2 = plt.subplot2grid((82, 116), (7, 65), colspan=51, rowspan=29)
            ax3 = plt.subplot2grid((82, 116), (46, 2), colspan=51, rowspan=29)
            ax4 = plt.subplot2grid((82, 116), (46, 65), colspan=51, rowspan=29)

            # Plot degradation curves
            ax1.semilogx(strain_lvl_or, gg0_cor[k], color='m',
                         marker='o', label='Kishida')
            ax1.set_ylim([0.0, 1.0])
            ax1.set_xlim([1e-5, 10.0])
            ax1.set_xlabel('Shear strain (%)')
            ax1.set_ylabel('G/G0')
            ax1.legend(loc='lower left')
            ax1.set_title('Shear modulus degradation curve')

            # Plot stress-strain curves
            ax2.plot(strain_lvl_or, shst_cor[k] * 0.001, color='m',
                     marker='o', label='Kishida')
            ax2.set_xlabel('Shear strain (%)')
            ax2.set_ylabel('Shear stress (kPa)')
            ax2.set_xlim([0.0, 5.0])
            ax2.set_ylim(bottom=0)
            ax2.legend(loc='lower right')
            ax2.set_title('Backbone curve')

            # Plot damping curves
            ax3.semilogx(strain_lvl_or, dam_cor[k] * 100, color='m',
                         marker='o', label='Kishida')
            ax3.set_xlim([1e-5, 10.0])
            ax3.set_ylim([0.0, 100.0])
            ax3.set_xlabel('Shear strain (%)')
            ax3.set_ylabel('Damping (%)')
            ax3.legend(loc='upper left')
            ax3.set_title('Damping curve')

            # Plot theta_tau curves
            ax4.plot([], [], color='m',
                     marker='o', label='Kishida')
            ax4.set_xlabel('Shear strain (%)')
            ax4.set_ylabel('Theta tau (-)')
            ax4.set_xlim([0.0, 5.0])
            ax4.set_ylim([0.0, 1.0])
            ax4.legend(loc='upper right')
            ax4.set_title('Theta-tau curve')

            # Add plot title
            flarge = round(1.50 * plt.rcParams['font.size'])
            fsmall = round(0.85 * plt.rcParams['font.size'])
            title = "Soil Layer " + str(k + 1) + " - " + data.inp['Soil Type'][k]
            plt.suptitle(title, fontsize=flarge, fontweight='bold',
                         horizontalalignment='left',
                         verticalalignment='top', x=0.02, y=0.97)

            plt.figtext(0.98, 0.97, txt, fontsize=fsmall,
                        horizontalalignment='right',
                        verticalalignment='top')

            plt.figtext(0.02, 0.03, note, fontstyle='oblique',
                        horizontalalignment='left',
                        verticalalignment='bottom')

            # Save figures
            name = os.path.join(out_dir, title)
            plt.savefig(name)
            plt.close()

    return strain_lvl_or, gg0_cor, shst_cor, dam_cor