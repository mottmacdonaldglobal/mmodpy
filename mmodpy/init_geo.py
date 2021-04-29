# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 13:49:24 2018

@author: kevin.stanton

This script executes initialization functions and serves as a constructor for 'geo.py'
"""


class update():

    def __init__(self, params=['Input.xlsx', '']):

        # Import libraries
        import xlrd
        import pandas as pd
        import numpy as np
        import os
        import sys

        # Import system paths to python functions
        self.cwd = os.getcwd()
        sys.path.insert(0, self.cwd)

        # Specify row indices defining input data sections
        r_aType = 0  # Analysis Type
        r_sParam = r_aType + 7  # Global Soil Parameters
        r_bCurve = r_sParam + 3  # Backbone Curve Parameters
        r_SRA = r_bCurve + 9  # SRA Parameters
        r_bedrock = r_SRA + 5 # Bedrock parameters
        r_aP = r_bedrock + 5  # Pile Analysis Parameters
        r_sBlock = r_aP + 6  # 3D Soil Block
        r_dLabel = r_sBlock + 11  # LS-DYNA Labels
        r_pbInfo = r_dLabel + 5  # Project and Building Information
        r_msLevs = r_pbInfo + 7  # MAT_HYSTERETIC_SOIL Strain Levels
        
        # Initialize parameters
        self.input_file = params[0]
        self.th_dir = params[1]

        # Define full input file directory (exn) and base name (input_file)
        if os.path.isfile(self.input_file):
            self.exn = self.input_file
            self.input_file = os.path.basename(self.input_file)
        elif os.path.isfile(os.path.join(self.cwd, self.input_file)):
            self.exn = os.path.join(self.cwd, self.input_file)
            
        # Define the main working directory (..\geopy within site-packages)
        self.mwd = os.path.dirname(os.path.realpath(__file__))
        
        # Define analysis type to be performed
        self.back_gen = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aType + 1, 6).value)
        self.sra_an = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aType + 2, 6).value)
        self.csp_an = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aType + 3, 6).value)
        self.so_bl = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aType + 4, 6).value)
        self.rr_PlaneStrain = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aType + 5, 6).value)

        # Define general soil properties
        self.inp = pd.read_excel(self.exn, sheet_name='Soil')

        # Define depth of groundwater table (m)
        self.depth_wt = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sParam + 1, 6).value)

        # Define Darendeli model to be used as target
        self.darend_soil = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_bCurve + 1, 6).value)
        if self.back_gen == 'Generated':
            # Define whether strain level for strength constrains are after
            # Vardanega
            self.str_vard = str(
                xlrd.open_workbook(
                    self.exn).sheet_by_name('Generic').cell(
                    r_bCurve + 2, 6).value)
            # Define whether strain level for strength constrains are after
            # Hayashi
            self.str_haya = str(
                xlrd.open_workbook(
                    self.exn).sheet_by_name('Generic').cell(
                    r_bCurve + 3, 6).value)
            # Define strain level for strength constrains (cohesive soils)
            self.str_match_clay = float(
                xlrd.open_workbook(
                    self.exn). sheet_by_name('Generic').cell(
                    r_bCurve + 4, 6).value)
            # Define strain level for strength constrains (cohesionless soils)
            self.str_match_sand = float(
                xlrd.open_workbook(
                    self.exn). sheet_by_name('Generic').cell(
                    r_bCurve + 5, 6).value)
            # Define strength ratio for strength constrains (cohesive soils)
            self.str_ratio_clay = float(
                xlrd.open_workbook(
                    self.exn). sheet_by_name('Generic').cell(
                    r_bCurve + 6, 6).value)
            # Define strength ratio for strength constrains (cohesionless
            # soils)
            self.str_ratio_sand = float(
                xlrd.open_workbook(
                    self.exn). sheet_by_name('Generic').cell(
                    r_bCurve + 7, 6).value)

        # SRA input parameters
        if self.sra_an == "Yes":
            # Define number of time histories considered
            if xlrd.open_workbook(self.exn).sheet_by_name('Generic').cell(
                    r_SRA + 1, 6).value == 'All (1D)':
                self.num_mot = int(len([name for name in os.listdir(self.th_dir) if os.path.isfile(os.path.join(self.th_dir, name))]))
                self.th_mode = 'All (1D)'
            elif xlrd.open_workbook(self.exn).sheet_by_name('Generic').cell(
                    r_SRA + 1, 6).value == 'All (2D)':
                self.num_mot = int(len([name for name in os.listdir(self.th_dir) if os.path.isfile(os.path.join(self.th_dir, name))]))
                self.th_mode = 'All (2D)'
            else:
                self.num_mot = int(xlrd.open_workbook(
                        self.exn).sheet_by_name('Generic').cell(
                                r_SRA + 1, 6).value)
                self.th_mode = 'default'
            if self.num_mot == 0:
                raise AssertionError('Directory for input time-histories contains no data.')
            ## Modify time history inputs if running 1D only
            # if self.th_mode == 'All (1D)':
            #     # sourcefiles = os.listdir(self.th_dir)
            #     destinationpath = os.path.join(self.mwd,'Time_Histories_1D_xyDuplicates')
            #     if os.path.isdir(destinationpath):
            #         shutil.rmtree(destinationpath)
            #     os.makedirs(destinationpath)
            #     i = 1
            #     for file in sourcefiles:
            #         newName_x = 'TH' + str(i) + '.csv'
            #         newName_y = 'TH' + str(i+1) + '.csv'
            #         if file.endswith('.csv'):
            #             shutil.copy(os.path.join(self.th_dir, file), os.path.join(destinationpath, newName_x))
            #             tempData = pd.read_csv((os.path.join(self.th_dir, file)),header=None)
            #             tempData.loc[0,0] = tempData.loc[0,0] + '_y'
            #             tempData.loc[0,1] = ''
            #             tempData.to_csv(os.path.join(destinationpath, newName_y),header=False,index=False)
            #         i = i + 2
            #     if os.path.isdir(os.path.join(self.th_dir,'Time_Histories_1D_xyDuplicates')):
            #         shutil.rmtree(os.path.join(self.th_dir,'Time_Histories_1D_xyDuplicates'))
            #     copy_tree(destinationpath, os.path.join(self.th_dir,'Time_Histories_1D_xyDuplicates'))
            #     self.th_dir = destinationpath
            #     self.num_mot = int(self.num_mot * 2)
            # Scale factor for time histories based on hazard map
            self.sc_factor = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_SRA + 2, 6).value
            # Viscous damping for SIREN SRA (default = 0.005)
            self.vis_damp = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_SRA + 3, 6).value
            self.logplot = False  # Flag to set linear or log scale
        else:
            self.vis_damp = 0.001

        # Bedrock input parameters
        if self.sra_an == 'Yes' or self.csp_an == 'Yes' or self.so_bl == 'Yes':
            # Bedrock shear wave velocity (m/sec)
            self.bed_sh = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_bedrock + 1, 6).value
            # Bedrock compression wave velocity (m/sec)
            self.bed_cp = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_bedrock + 2, 6).value
            # Bedrock density (kg/m3)
            self.bed_de = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_bedrock + 3, 6).value

        # Define ALP and LS-DYNA pile analysis parameters
        if self.csp_an == "Yes":
            # Read the number of piles to be analysed
            self.num_piles = int(
                max(pd.read_excel(self.exn, sheet_name='Pile', header=0).iloc[:, 0]))
            # Define the pile properties from the input data
            self.pile_prop = pd.read_excel(
                self.exn, sheet_name='Pile', header=0).iloc[0:self.num_piles, :]
            # Maximum beam element length in ALP (m)
            self.alp_len = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aP + 1, 6).value
            # Maximum beam element length in DYNA (m)
            self.dyna_len = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aP + 2, 6).value
            # Size of the finite element verification model (m)
            self.model_size = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aP + 3, 6).value
            # Size of the solid elements in verification model (m)
            self.mesh_size = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_aP + 4, 6).value

        # Define LS-DYNA soil block dimensions and element size
        if self.so_bl == "Yes":
            # Read the number of shallow foundations to be analysed
            self.num_foot = int(max(pd.read_excel(
                self.exn, sheet_name='Footings and Grade Beams', header=0).iloc[:, 0]))
            # Define the shallow foundation properties from the input data
            self.foot_prop = pd.read_excel(
                self.exn, sheet_name='Footings and Grade Beams', header=0).iloc[0:self.num_foot, :]
            # Outer dimensions of foundations in global x-direction (m)
            self.lx_found = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 1, 6).value
            # Outer dimensions of foundations in global y-direction (m)
            self.ly_found = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 2, 6).value
            # Scale factor for damping value in x - direction (-)
            self.dxfac = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 3, 6).value
            # Scale factor for damping value in y - direction (-)
            self.dyfac = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 4, 6).value
            # Depth of bottom surface grade beam (m)
            self.excv_dep = str(
                xlrd.open_workbook(
                    self.exn).sheet_by_name('Generic').cell(
                    r_sBlock + 5, 6).value)
            # Extent of soil block in global x-direction (m)
            self.lx_bl = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 6, 6).value
            # Extent of soil block in global y-direction (m)
            self.ly_bl = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 7, 6).value
            # Maximum element size in vicinity of foundations (m)
            self.elsiz_bl = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 8, 6).value
            # Maximum element size at far-field boundary (m)
            self.max_elsiz_bl = xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_sBlock + 9, 6).value

        # Define LS-DYNA label ranges
        # Start of soil include label range (sra.key, soil.key)
        self.s_id = int(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_dLabel + 1, 6).value)
        # Start of pile include label range (csp.key)
        self.p_id = int(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_dLabel + 2, 6).value)
        # Increment for pile include label range (csp#.key)
        self.p_id_incr = int(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_dLabel + 3, 6).value)

        # Define project information
        # Define job number
        self.jn = str(int(xlrd.open_workbook(self.exn).sheet_by_name(
            'Generic').cell(r_pbInfo + 1, 6).value))
        # Define job title
        self.jt = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_pbInfo + 2, 6).value)
        # Define job location
        self.jl = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_pbInfo + 3, 6).value)
        # Define job subtitle
        self.js = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_pbInfo + 4, 6).value)
        # Define soil properties
        self.sp = str(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_pbInfo + 5, 6).value)

        # Define strain levels for LS-DYNA analysis
        str1 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 1,
                6).value)
        str2 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 2,
                6).value)
        str3 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 3,
                6).value)
        str4 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 4,
                6).value)
        str5 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 5,
                6).value)
        str6 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 6,
                6).value)
        str7 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 7,
                6).value)
        str8 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 8,
                6).value)
        str9 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 9,
                6).value)
        str10 = float(
            xlrd.open_workbook(
                self.exn).sheet_by_name('Generic').cell(
                r_msLevs + 10,
                6).value)
        self.strain_lvl_dy = np.array(
            [str1, str2, str3, str4, str5, str6, str7, str8, str9, str10])
