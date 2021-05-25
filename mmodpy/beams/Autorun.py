# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 12:39:44 2018

@author: kevin.stanton

This class uses the autopy function to run an XTRACT file and export the
analysis results. The output file is a ".txt" file that is saved and moved in
the current working directory.

Notes:
    Developed on a 1920x1080 display. All mouse movements are normalized and
    scaled by these values.
"""

import tkinter as tk
import subprocess
from autopy3 import key
import time
import psutil
import os
from .util import util


class Autorun:

    def __init__(self, aData, mcSettings, xpj_name=''):

        # Get current working directory and screen size
        self.cwd = str(aData.pSave)

        # Get screen dimensions
        root = tk.Tk()
        self.sWidth = root.winfo_screenwidth()
        self.sHeight = root.winfo_screenheight()
        self.sWidthF = self.sWidth / 1920
        self.sHeightF = self.sHeight / 1080
        root.quit()

        # Get settings data
        self.mcSettings = mcSettings

        # Get autodyna module path
        import autodyna.Resources.XTRACT as pathIndex
        self.xtractPath = os.path.dirname(pathIndex.__file__)
        
        # Determine path to .xpj file to run
        if xpj_name == '':
            self.xpj_name = self.cwd + r'\\' + str(aData.file_name) + '.xpj'
        else:
            self.xpj_name = os.path.join(os.getcwd(), xpj_name)

    def __launch_app__(self, aData):

        # Launch the application
        try:
            p = psutil.Popen([aData.pXTRACT, self.xpj_name],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        except FileNotFoundError:
            try:
                p = psutil.Popen(['C:\\Program Files (x86)\\TRC\\XTRACT\\XTRACT.exe',
                                  self.xpj_name],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            except FileNotFoundError:
                try:
                    p = psutil.Popen(['C:\\ProgramData\\XTRACT\\XTRACT.exe',
                                      self.xpj_name],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                except FileNotFoundError:
                    try:
                        p = psutil.Popen([os.path.join(self.xtractPath,
                                                       'XTRACT.exe'),
                                          self.xpj_name],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
                    except FileNotFoundError:
                        raise AssertionError(
                            'Error: XTRACT not properly installed/associated with wrong path')

        # Assign subproccess
        self.proc = p
        util.confirm(aData)
        if not aData.mInter:
            time.sleep(self.mcSettings.sleepLaunch)

    def __run_anal__(self, aData, n_strainlim):
        time.sleep(2)
        util.click(1, 1)
        time.sleep(1)
        util.click(1, 1)
        time.sleep(1)
        key.tap(key.K_F5)  # Run analysis
        print('\n   Running XTRACT...')
        util.confirm(aData)

        if aData.timingControl == 'automatic':
            # CPU usage tracking
            time.sleep(15)
            cpu_hist, stop_flag = util.sample_cpu_percent(self.proc)
            time.sleep(2)
        else:
            # Fixed analysis timing
            c = 1  # Reduce run time if plots are not shown
            if aData.showPlots:
                c = 0.9
            analysis_time = c * (self.mcSettings.sleepAt1 * 2 + float(n_strainlim) * self.mcSettings.sleepAt2) * \
                float(aData.section_area * self.mcSettings.sleepAt3 * 1e6) / float(600 * 600)
            if not aData.mInter:
                print('\n            sleeping for ' +
                      str(round(analysis_time +
                                self.mcSettings.sleepRunAnal, 1)) +
                      ' seconds \n')
                time.sleep(analysis_time)
                time.sleep(self.mcSettings.sleepRunAnal)

    def __show_res__(self, aData, n_strainlim):

        util.click(1, 1)  # Upper-left corner to get to window
        time.sleep(1)
        util.click(1, 1)
        time.sleep(1)
        key.toggle('d', True, key.MOD_CONTROL)
        key.toggle('d', False, key.MOD_CONTROL)
        time.sleep(self.mcSettings.sleepSLc1 + self.mcSettings.sleepSLc2)
        util.click(int(self.mcSettings.clickPosition1 *
                       self.sWidthF), int(self.mcSettings.clickPosition2 *
                                          self.sHeightF), 'right') # Click 'Project Name'
        time.sleep(self.mcSettings.sleepSLc1 + self.mcSettings.sleepSLc2)
        util.click(int(self.mcSettings.clickPosition3 *
                       self.sWidthF), int(self.mcSettings.clickPosition4 *
                                          self.sHeightF)) # Click 'Show Data'
        util.confirm(aData)
        if not aData.mInter:
            time.sleep(self.mcSettings.sleepSLc1 + self.mcSettings.sleepSLc2)

    def __select_data__(self, aData, n_strainlim):

        # Add cwd to clipboard
        saveName = os.path.join(self.cwd,str(aData.file_name) + '.txt')
        util.addToClipBoard(saveName)

        # Select and save data
        util.click(int(self.mcSettings.clickPosition5 *
                       self.sWidthF), int(self.mcSettings.clickPosition6 *
                                          self.sHeightF)) # Click 'Select All'
        time.sleep(self.mcSettings.sleepSLc1 + self.mcSettings.sleepSLc2)
        util.click(int(self.mcSettings.clickPosition7 *
                       self.sWidthF), int(self.mcSettings.clickPosition6 *
                                          self.sHeightF)) # Click 'Save Data'
        time.sleep(self.mcSettings.sleepSLc1 + self.mcSettings.sleepSLc2 * 0.5)
        key.toggle('v', True, key.MOD_CONTROL)  # Paste cwd string
        key.toggle('v', False, key.MOD_CONTROL)
        time.sleep(self.mcSettings.sleepPaste)
        key.tap(key.K_RIGHT) # Move to end of pasted string
        if self.mcSettings.controlMode == 'default':
            key.tap(key.K_BACKSPACE) # Clear trailing space
        key.tap(key.K_RETURN) # Execute cd command
        time.sleep(self.mcSettings.sleepSLc1 + self.mcSettings.sleepSLc2)
        util.click(int(self.mcSettings.clickPosition8 *
                       self.sWidthF), int(self.mcSettings.clickPosition9 *
                                          self.sHeightF)) # Select dialog box at outermost level
        key.tap(key.K_RETURN)  # Save
        time.sleep(self.mcSettings.sleepSave)

    def __close_app__(self):

        # Close the program
        os.system("TASKKILL /F /IM xtract.exe")
        time.sleep(1)

    def autorun_xtract(self, aData, n_strainlim=1):

        # Autorun XTRACT analysis
        self.__launch_app__(aData)
        self.__run_anal__(aData, n_strainlim)
        self.__show_res__(aData, n_strainlim)
        self.__select_data__(aData, n_strainlim)
        saveName = os.path.join(self.cwd,str(aData.file_name) + '.txt')

        # Check for errors and correct if possible
        if os.path.isfile(saveName) == False:
            time.sleep(10)
            if os.path.isfile(saveName):
                self.__show_res__(aData, n_strainlim)
                self.__select_data__(aData, n_strainlim)
            else:
                time.sleep(20)
                if os.path.isfile(saveName):
                    self.__show_res__(aData, n_strainlim)
                    self.__select_data__(aData, n_strainlim)
                else:
                    time.sleep(20)
                    self.__show_res__(aData, n_strainlim)
                    self.__select_data__(aData, n_strainlim)

        self.__close_app__()
