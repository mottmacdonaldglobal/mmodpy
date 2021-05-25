# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 11:52:55 2018

@author: kevin.stanton

Contains various mouse control and timing settings for hijacking XTRACT.
"""

import pandas as pd
import os

class settings:

    def __init__(self, aData):

        # Mouse position settings (change if mouse positions are incorrect)
        self.clickPosition1 = 60 # Click 'Project Name' - x
        self.clickPosition2 = 105 # Click 'Project Name' - y
        self.clickPosition3 = 125 # Click 'Show Data' - x
        self.clickPosition4 = 170 # Click 'Show Data' - y
        self.clickPosition5 = 700 # Click 'Select All' - x
        self.clickPosition6 = 1025 # Click 'Select All' or 'Save Data' - y
        self.clickPosition7 = 842 # Click 'Save Data' - x
        self.clickPosition8 = 1300 # Select dialog box at outermost level - x
        self.clickPosition9 = 400 # Select dialog box at outermost level - y
        
        # Time settings (all wait times in seconds, changes not necessary)
        self.sleepLaunch = 6  # Post-launch wait time
        self.sleepRunAnal = 4  # Post-execute XTRACT analysis command wait time
        self.sleepAt = 1.2  # Proxy for dynamically assigned additional sleep time while analysis is running
        self.sleepSLc = 2  # Proxy for data extraction sleep time
        self.sleepPaste = 1.0  # Post-paste cwd in Windows Explorer wait time
        self.sleepSave = 1.5  # Post-execute save command wait time

        # Compute timing constants (changes not necessary)
        self.sleepAt1 = 7 * self.sleepAt  # Analysis timing constant 1
        self.sleepAt2 = 5 * self.sleepAt  # Analysis timing constant 2
        self.sleepAt3 = 1.6 * self.sleepAt  # Analysis timing constant 3
        self.sleepSLc1 = 1.5 * self.sleepSLc  # Data extraction contant 1
        self.sleepSLc2 = 0.1 * self.sleepSLc  # Data extraction contant 2
        
        # Read in copy/paste behavior handling option
        self.controlMode = aData.controlMode 
        
        # Read in user-defined mouse click positions if specified
        if os.path.isfile(aData.settings_path):
            clickPositions = pd.read_csv(aData.settings_path, dtype=str, 
                                              header=None, sep=',\s+', 
                                              engine='python', index_col=0,
                                              comment='#')
            self.clickPosition1 = int(clickPositions.loc['clickPosition1']) # Click 'Project Name' - x
            self.clickPosition2 = int(clickPositions.loc['clickPosition2']) # Click 'Project Name' - y
            self.clickPosition3 = int(clickPositions.loc['clickPosition3']) # Click 'Show Data' - x
            self.clickPosition4 = int(clickPositions.loc['clickPosition4']) # Click 'Show Data' - y
            self.clickPosition5 = int(clickPositions.loc['clickPosition5']) # Click 'Select All' - x
            self.clickPosition6 = int(clickPositions.loc['clickPosition6']) # Click 'Select All' or 'Save Data' - y
            self.clickPosition7 = int(clickPositions.loc['clickPosition7']) # Click 'Save Data' - x
            self.clickPosition8 = int(clickPositions.loc['clickPosition8']) # Select dialog box at outermost level - x
            self.clickPosition9 = int(clickPositions.loc['clickPosition9']) # Select dialog box at outermost level - y
            

        