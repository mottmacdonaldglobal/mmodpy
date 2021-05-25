"""
Created on Thu Apr 26 13:01:57 2018

@author: kevin.stanton

This class provides helper functions for autodyna_ac
"""

import pandas as pd
import numpy as np
import math
import time
from scipy.interpolate import interp1d
import os
from tkinter import messagebox
import win32api
import win32con


class util:

    def confirm(aData):

        if aData.mInter:
            answer = messagebox.askyesno(
                "Breakpoint Reached", "Ready to continue?")
            if not answer:
                raise AssertionError("Process terminated by user.")

    def click(x, y, button='left'):

        win32api.SetCursorPos((x, y))

        if button == 'right':
            win32api.mouse_event(win32con.MOUSEEVENTF_RIGHTDOWN, x, y, 0, 0)
            win32api.mouse_event(win32con.MOUSEEVENTF_RIGHTUP, x, y, 0, 0)
        elif button == 'left':
            win32api.mouse_event(win32con.MOUSEEVENTF_LEFTDOWN, x, y, 0, 0)
            win32api.mouse_event(win32con.MOUSEEVENTF_LEFTUP, x, y, 0, 0)

    def move(x, y):

        win32api.SetCursorPos((x, y))

    def combineKeys(keyNames, outFile='autoBeams.key'):
        """
        This function combines .key files present in the current woking directory
        into one file named 'autoBeams.key'.

        Parameters
        ----------
            keyNames : list
                - file names or full file paths to be combined as a list of strings
                - list elements must represent existing .key files
                - example: ['file1.key', 'file2.key', ...]
            outFile : str
                - optional [default='autoBeams.key']
                - file name to write combined .key data

        Returns
        -------
            [outfile].key : file
                - .key with data combined from keyNames list
                - written to the current working directory
        """

        # Define list of terms to delete in order to facilitate the merge
        deleteList = ['*KEYWORD', '*END']
        # Write the raw combined file (temp)
        with open('temp.key', 'wb') as outfile:
            for i in keyNames:
                with open(i, 'rb') as infile:
                    outfile.write(infile.read())
        # Remove the required strings and write the final file
        with open('temp.key') as oldFile, open(outFile, 'w') as newFile:
            newFile.write('*KEYWORD\n')
            for line in oldFile:
                if not any(
                        badString in line for badString in deleteList) and line.strip():
                    newFile.write(line)
            newFile.write('*END')
        # Delete the temp file
        os.remove('temp.key')

    def cutOutlier(x, y):

        numCut = 0
        i = 0
        xMod = []
        yMod = []
        while i < len(y) - 1:
            if float(y[i]) < float(y[i + 1]):
                xMod.append(x[i])
                yMod.append(y[i])
            elif i > 9:
                xMod.append(x[i])
                yMod.append(y[i])
            else:
                numCut = numCut + 1
            i = i + 1

        return xMod, yMod, numCut

    def expYield(M_phi, aData):

        steelThreshold = float(aData.Fy) / float(aData.Es)
        print('         Concrete yield strain = ' + str(aData.ey))
        print('         Steel yield strain = ' + str(steelThreshold))

        # Read in data
        Mxx_raw, phi_xx_raw, cMin_xx, sMin_xx, sMax_xx = M_phi.loc[:, 'Mxx (N-m)'], \
            M_phi.loc[:, 'phi_xx (1/m)'], M_phi.loc[:, 'ConfinedMinStrain_xx'], \
            M_phi.loc[:, 'SteelMinStrain_xx'], M_phi.loc[:, 'SteelMaxStrain_xx']

        Myy_raw, phi_yy_raw, cMin_yy, sMin_yy, sMax_yy = M_phi.loc[:, 'Myy (N-m)'], \
            M_phi.loc[:, 'phi_yy (1/m)'], M_phi.loc[:, 'ConfinedMinStrain_yy'], \
            M_phi.loc[:, 'SteelMinStrain_yy'], M_phi.loc[:, 'SteelMaxStrain_yy']

        # Remove outliers
        phi_xx, Mxx, numCutX = util.cutOutlier(phi_xx_raw, Mxx_raw)
        phi_yy, Myy, numCutY = util.cutOutlier(phi_yy_raw, Myy_raw)

        # Find coordinate of My
        # xx direction
        i = 1
        while i <= len(Mxx_raw):
            index = 0
            if abs(float(cMin_xx[i])) > float(aData.ey):
                print('         Governing yield - xx = confined concrete (XTRACT min)')
                index = min(max(i - 1, 4), 10)
                i = len(Mxx_raw)
            elif abs(float(sMax_xx[i])) > float(aData.Fy) / float(aData.Es):
                print('         Governing yield - xx = steel (XTRACT max)')
                index = min(max(i - 1, 4), 10)
                i = len(Mxx_raw)
            elif abs(float(sMin_xx[i])) > float(aData.Fy) / float(aData.Es):
                print('         Governing yield - xx = steel (XTRACT min)')
                index = min(max(i - 1, 4), 10)
                i = len(Mxx_raw)
            i = i + 1

        # Check that an appropriate index was found
        if index == 0:
            raise AssertionError(
                'Yield moment cannot be interpreted! Moment curvature must be computed for larger plastic rotations.')
        print('         index_xx =' + str(index))

        ind_xx = index
        xP_xx = float(phi_xx_raw[index])
        yP_xx = float(Mxx_raw[index])

        # yy direction
        i = 1
        while i <= len(Myy_raw):
            index = 0
            if abs(float(cMin_yy[i])) > float(aData.ey):
                print('         Governing yield - yy = confined concrete (XTRACT min)')
                index = min(max(i - 1, 4), 10)
                i = len(Myy_raw)
            elif abs(float(sMax_yy[i])) > float(aData.Fy) / float(aData.Es):
                print('         Governing yield - yy = steel (XTRACT max)')
                index = min(max(i - 1, 4), 10)
                i = len(Myy_raw)
            elif abs(float(sMin_yy[i])) > float(aData.Fy) / float(aData.Es):
                print('         Governing yield - yy = steel (XTRACT min)')
                index = min(max(i - 1, 4), 10)
                i = len(Myy_raw)
            i = i + 1

        # Check that an appropriate index was found
        if index == 0:
            raise AssertionError(
                'Yield moment cannot be interpreted! Moment curvature must be computed for larger plastic rotations.')
        print('         index_yy =' + str(index))

        # Determine outputs for the computed index
        ind_yy = index
        xP_yy = float(phi_yy_raw[index])
        yP_yy = float(Myy_raw[index])

        # Determine additional outputs
        M_ult_xx = float(max(Mxx))
        M_ult_yy = float(max(Myy))
        k_ult_xx = np.interp(M_ult_xx, Mxx, phi_xx)
        k_ult_yy = np.interp(M_ult_yy, Myy, phi_yy)

        # Combine results for x and y
        xP = np.array([abs(xP_xx), abs(xP_yy)])  # phi
        yP = np.array([abs(yP_xx), abs(yP_yy)])  # M
        M_ult = np.array([abs(M_ult_xx), abs(M_ult_yy)])
        k_ult = np.array([abs(k_ult_xx), abs(k_ult_yy)])
        ind = np.array([abs(ind_xx), abs(ind_yy)])

        return xP, yP, M_ult, k_ult, ind

    def getName(aData, filetype='title'):
        partType = 'Beam-RC'
        secType = aData.section_type
        if secType == "Rectangular Beam" or secType == 'T Beam':
            secType = 'Beam'
        elif secType == 'Circular' or secType == 'Rectangular Column':
            secType = 'Column'
        if filetype == '.key':  # .key file name
            filename = str(partType + '_' + secType + '_' +
                           str(aData.file_name) + '.key')
        else:  # PART naming within PRIMER
            filename = str(partType + '_' + secType +
                           ';' + str(aData.file_name))
        return filename

    def getPM(aData, n_strainlim, npoints):
        filename = aData.file_name
        pmxx_array, pmyy_array = pd.DataFrame(), pd.DataFrame()
        dat = pd.read_csv(aData.pSave + '\\' + str(filename) + '.txt',
                          delim_whitespace=True,
                          error_bad_lines=False,
                          skiprows=15)
        colstep = 38  # 38 columns added for each new PM curve
        colind = [12, 14, 31, 34]
        for i in range(0, n_strainlim):
            newind = [x + (colstep * i) for x in colind]
            pmxx, pmyy = dat.iloc[:, newind[0:2]], dat.iloc[:, newind[2:4]]
            pmxx.columns, pmyy.columns = [
                'P' + str(i), 'M' + str(i)], ['P' + str(i), 'M' + str(i)]
            pmxx_array = pd.concat([pmxx_array, pmxx], axis=1)
            pmyy_array = pd.concat([pmyy_array, pmyy], axis=1)
        pmxx_fin = pd.DataFrame(data=None, columns=['Px', 'Mx'])
        pmyy_fin = pd.DataFrame(data=None, columns=['Py', 'My'])

        # get envelope
        pmxx_fin = util.getEnvelope(pmxx_array, npoints, axis='xx')
        pmyy_fin = util.getEnvelope(pmyy_array, npoints, axis='yy')
        pm_envelope = pd.concat([pmxx_fin, pmyy_fin], axis=1)
        pm_all = pd.concat([pmxx_array, pmyy_array], axis=1)
        pm_all.to_csv(
            aData.pSave +
            '\\' +
            'Moment_Interaction' +
            filename +
            '.csv',
            index=False)
        pm_envelope.to_csv(aData.pSave + '\\' + 'Moment_Interaction_Envelope' +
                           filename + '.csv', index=False)
        pmxx_array.to_csv(
            aData.pSave +
            '\\' +
            'PMxx' +
            filename +
            '.csv',
            index=False)
        pmyy_array.to_csv(
            aData.pSave +
            '\\' +
            'PMyy' +
            filename +
            '.csv',
            index=False)

        return pm_envelope

    def getPlasticHingeLength(aData):
        shape = aData.section_type
        # Check plastic hinge length factor
        phl_xx = float(aData.phl_xx)
        phl_yy = float(aData.phl_yy)
        if float(phl_xx) < -1 or float(phl_yy) < -1:
            print('========================================')
            print('WARNING: plastic hinge length factor > 1')
            print('========================================')

        if shape == 'Circular':
            D = float(aData.diameter_circ) / 1000.0
            if phl_xx == 0 or phl_xx < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_xx = D / 2  # default equations
            elif phl_xx > 0:  # absolute (return in m)
                pl_hi_xx = phl_xx / 1000.0
            elif phl_xx < 0:  # factor
                pl_hi_xx = D * abs(phl_xx)

            if phl_yy == 0 or phl_yy < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_yy = D / 2  # defaulte quations
            elif phl_yy > 0:  # absolute
                pl_hi_yy = phl_yy / 1000.0
            elif phl_yy < 0:  # factor
                pl_hi_yy = D * abs(phl_yy)

        elif (shape == 'Rectangular Beam') | (shape == 'Rectangular Column'):
            # Get the width of the section (m)
            w = float(aData.rect_width) / 1000
            # Get the height of the section (m)
            h = float(aData.rect_height) / 1000
            # Plastic hinge length (m)
            if phl_xx == 0 or phl_xx < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_xx = h / 2  # default equations
            elif phl_xx > 0:  # absolute (return in m)
                pl_hi_xx = phl_xx / 1000.0
            elif phl_xx < 0:  # factor
                pl_hi_xx = h * abs(phl_xx)

            if phl_yy == 0 or phl_yy < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_yy = w / 2
            elif phl_yy > 0:  # absolute (return in m)
                pl_hi_yy = phl_yy / 1000.0
            elif phl_yy < 0:  # factor
                pl_hi_yy = w * abs(phl_yy)

        elif shape == 'T Beam':
            # Get the width of the flange (m)
            B = float(aData.f_width) / 1000
            # Get the height of the flange (m)
            h = float(aData.f_thick) / 1000
            # Get the width of the web (m)
            b = float(aData.t_width) / 1000
            # Get the height of the web (m)
            H = float(aData.t_height) / 1000 - float(aData.f_thick) / 1000
            if phl_xx == 0 or phl_xx < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_xx = (H + h) / 2  # default equations
            elif phl_xx > 0:  # absolute (return in m)
                pl_hi_xx = phl_xx / 1000.0
            elif phl_xx < 0:  # factor
                pl_hi_xx = (H + h) * abs(phl_xx)

            if phl_yy == 0 or phl_yy < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_yy = b / 2  # default equations
            elif phl_yy > 0:  # absolute (return in m)
                pl_hi_yy = phl_yy / 1000.0
            elif phl_yy < 0:  # factor
                pl_hi_yy = (b) * abs(phl_yy)

        elif shape == 'Cruciform':
            # Total width (m)
            bX = float(aData.bX) / 1000
            # Total height (m)
            bY = float(aData.bY) / 1000
            if phl_xx == 0 or phl_xx < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_xx = bX / 2  # default equations
            elif phl_xx > 0:  # absolute (return in m)
                pl_hi_xx = phl_xx / 1000.0
            elif phl_xx < 0:  # factor
                pl_hi_xx = bX * abs(phl_xx)

            if phl_yy == 0 or phl_yy < -1:  # <-1 is a safety measure for factors > 1
                pl_hi_yy = bY / 2  # default equations
            elif phl_yy > 0:  # absolute (return in m)
                pl_hi_yy = phl_yy / 1000.0
            elif phl_yy < 0:  # factor
                pl_hi_yy = bY * abs(phl_yy)

        return pl_hi_xx, pl_hi_yy

    def getNormalize(dat):
        # get normalization factors
        pmax, pmin, mmax, mmin = 0, 0, 0, 0
        for i in range(0, len(dat.columns), 2):
            pmax_i = max(list(dat[dat.columns[i]]))
            mmax_i = max(list(dat[dat.columns[i + 1]]))
            pmin_i = min(list(dat[dat.columns[i]]))
            mmin_i = min(list(dat[dat.columns[i + 1]]))
            if pmax_i > pmax:
                pmax = pmax_i
            if mmax_i > mmax:
                mmax = mmax_i
            if pmin_i < pmin:
                pmin = pmin_i
            if mmin_i < mmin:
                mmin = mmin_i
        plen = abs(pmax) + abs(pmin)
        mlen = abs(mmax) + abs(mmin)

        return plen, mlen

    def getEnvelope(dat, nsteps, axis):
        # dat is the pandas dataframe that contains
        # the combined PM interaction curves for all limiting strains tested

        # nsteps is how many points you want on the final PM curve

        # step 1 - find max values to get the normalization values
        plen, mlen = util.getNormalize(dat)

        # step 2 - initialize the angles at which you want to generate the PM
        # envelope
        rad_vec = np.linspace(-np.pi, np.pi, nsteps)
        polar_array = np.zeros(shape=(nsteps, int(len(dat.columns) / 2.0 + 1)))
        polar_array[:, 0] = rad_vec

        # step 3 - loop through each limiting strain curve
        for i in range(0, len(dat.columns), 2):
            # normalize each limiting strain curve
            p = np.divide(list(dat[(dat.columns[i])]), plen)
            m = np.divide(list(dat[(dat.columns[i + 1])]), mlen)

            # polar coordinates
            rad = [math.atan2(p[x], m[x]) for x in range(0, len(p))]
            mag = [math.sqrt(pow(p[x], 2) + pow(m[x], 2))
                   for x in range(0, len(p))]
            polar_vec = pd.DataFrame({'T': rad, 'R': mag}).reindex(
                columns=['T', 'R'])  # constrain order of columns
            polar_sort = np.array(polar_vec.sort_values(by=['T']))
            polar_interp = interp1d(polar_sort[:,
                                               0],
                                    polar_sort[:,
                                               1],
                                    bounds_error=False,
                                    fill_value="extrapolate")(rad_vec)
            polar_array[:, int(i / 2.0 + 1)] = polar_interp

        # step 4 - now that we have polar coords for all curves at the same angles decided by rad_vec
        # compare the magnitude at each angle then pick max

        finalpm = np.zeros(shape=(nsteps, 2))  # 2 columns for P, M
        for r in range(0, len(polar_array[:, 0])):
            row = list(polar_array[r, 1:])
            # index of curve corresponding to the envelope curve
            ind = row.index(max(row))
            # now re-calculate the actual, non-normalized P and M values from
            # normalized polar coords
            p = plen * polar_array[r, ind + 1] * math.sin(polar_array[r, 0])
            m = mlen * polar_array[r, ind + 1] * math.cos(polar_array[r, 0])
            finalpm[r, :] = [m, p]
        if axis == 'xx':
            finalpm_df = pd.DataFrame(
                {'Px': finalpm[:, 1], 'Mx': finalpm[:, 0]})
        else:  # yy
            finalpm_df = pd.DataFrame(
                {'Py': finalpm[:, 1], 'My': finalpm[:, 0]})

        return finalpm_df

    def cpu_percent(proc):
        """CPU check for timing analysis"""

        # Check CPU usage of a process
        cpu_pct = 0.0
        cpu_pct = proc.cpu_percent()

        return cpu_pct

    def sample_cpu_percent(proc):
        """Flag when analysis is finished"""

        # Initialize
        stop_flag = 0
        cpu_hist = 0.0
        j = 0

        while j < 10000:
            cpu_pct = util.cpu_percent(proc)
            time.sleep(1)
            cpu_hist = max(cpu_pct, cpu_hist)
            j += 1

            if j > 10:
                if cpu_pct < 3:
                    stop_flag = 1
                    print('            Xtract analysis complete \n')
                    break

        return cpu_hist, stop_flag

    def addToClipBoard(text):

        df = pd.DataFrame([text])
        df.to_clipboard(index=False, header=False)
