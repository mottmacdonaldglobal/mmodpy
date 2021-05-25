class info_hysBeams():

    def __init__(self, data, iRun):

        ##################################################
        ### Static data (same for all batch instances) ###
        ##################################################

        # General project information
        self.author = data[1][0]  # Name of author
        self.company = data[1][1]  # Company name
        self.job_name = data[1][2]  # Project name

        # Specify XTRACT path
        # Default = 'C:\Program Files (x86)\TRC\XTRACT\XTRACT.exe'
        self.pXTRACT = data[1][3]

        # Manual intervention mode
        self.mInter = data[1][4]

        # Batch mode options (False -> results only returned for given iSet)
        self.batchMode = data[1][5]
        
        # Input set for a single analysis (ignored if batchMode = True)
        self.iSet = data[1][6]
        self.iRun = int(iRun)
        
        # Option to write failure parameters
        self.fParams = data[1][7]
        
        # Treatment of T-beams
        self.tbOption = data[1][8]

        # Developer mode settings
        self.devMode = data[1][9]

        # Constraints on fitting for My interpretation ('strict' or 'flexible')
        self.cLimits = data[1][10]
        self.cleanCurves = data[1][11]

        # Flag to show plots in XTRACT analysis runs
        self.showPlots = data[1][12]

        # Read in control mode settings
        self.controlMode = data[1][13]

        # Code option
        self.codeOption = data[1][14]

        # Concrete factor
        self.concreteFactor = data[1][15]

        # Timing settings
        self.timingControl = data[1][16]

        # Specify path to save
        self.pSave = data[1][17]

        # Specify option to write terminal to log gile
        self.write_log = data[1][18]
        
        # User-defined mouse click positions file path
        self.settings_path = data[1][19]

        ########################################
        ### Generate Unique Names for XTRACT ###
        ########################################

        # CAUTION: Modify standard values below with care
        # In almost all cases, editing is not required
        self.mat_name_unc = 'AutoUnconfined'  # Unconfinced concrete
        self.mat_name_conc = 'AutoConfined'  # Confined concrete
        self.mat_name_steel = 'AutoSteel'  # Steel
        self.section_name = 'AutoSection'  # Section name
        self.material_core = 'AutoConfined'  # Concrete core material name
        self.material_cover = 'AutoUnconfined'  # Concrete cover material name
        self.material_steel = 'AutoSteel'  # Steel material name

        ###########################################################
        ### Array data (value required for each batch instance) ###
        ###########################################################

        # File name
        self.file_name = data[0].loc['info_filename', iRun]

        # Section, material and curve nomenclature
        try:
            # Name sections and curves based on Beam ID
            self.beam_num = int(self.file_name[2:])

        except (IndexError, ValueError) as error:
            # Name sections and curves based on run number
            self.beam_num = int(self.iRun + 1) * 10

        # Concrete material properties
        self.Fc = data[0].loc['Fc', iRun]  # 28-day compressive strength (MPa)
        self.Ft = data[0].loc['Ft', iRun]  # Tensile strength (MPa)
        self.ey = data[0].loc['ey', iRun]  # Yield strain
        self.ecu = data[0].loc['ecu', iRun]  # Crushing strain
        self.esp = data[0].loc['esp', iRun]  # Spalling strain
        self.ef = data[0].loc['ef', iRun]  # Failure strain
        self.Ec = data[0].loc['Ec', iRun]  # Modulus of elasticity (MPa)
        self.den = data[0].loc['den', iRun]  # Density of concrete (kg/m3)
        self.v = data[0].loc['v', iRun]  # Poisson's ratio of concrete

        # Steel material properties
        self.Fy = data[0].loc['Fy', iRun]  # Yield stress (MPa)
        self.Fu = data[0].loc['Fu', iRun]  # Fracture stress (MPa)
        self.esh = data[0].loc['esh', iRun]  # Strain at strain hardening
        self.esu = data[0].loc['esu', iRun]  # Failure strain
        self.Es = data[0].loc['Es', iRun]  # Modulus of elasticity (MPa)
        self.esy = float(self.Fy) / float(self.Es)  # Yield strain

        # General section properties
        # 'Rectangular Beam', 'Rectangular Column', 'Circular', 'T Beam', 'Cruciform'
        self.section_type = data[0].loc['section_type', iRun]  # Section type
        self.cover = data[0].loc['cover', iRun]  # concrete cover (mm)

        # Geometry for circular sections
        # Section Diameter (mm)
        self.diameter_circ = data[0].loc['diameter_circ', iRun]

        # Reinforcement properties for circular sections
        # Diameter of longitudinal bars (mm)
        self.bar_diam_circ = data[0].loc['bar_diam_circ', iRun]
        # Number of longitudinal bars
        self.num_bars_circ = data[0].loc['num_bars_circ', iRun]

        # Geometry for all rectangular sections
        self.rect_width = data[0].loc['rect_width', iRun]  # Section width (mm)
        # Section height (mm)
        self.rect_height = data[0].loc['rect_height', iRun]

        # Reinforcement properties for rectangular beam sections
        self.n_t = data[0].loc['n_t', iRun]  # Number of top bars
        self.n_b = data[0].loc['n_b', iRun]  # Number of bottom bars
        self.db_t = data[0].loc['db_t', iRun]  # Diameter of top bars (mm)
        self.db_b = data[0].loc['db_b', iRun]  # Diameter of bottom bars (mm)

        # Reinforcement properties for rectangular column sections
        # Number of bars ('4', '6', or '8')
        self.num_bars_col = data[0].loc['num_bars_col', iRun]
        # Diameter of reinforcement (mm)
        self.bar_diam_col = data[0].loc['bar_diam_col', iRun]

        # Geometry for T Beam sections
        self.t_width = data[0].loc['t_width', iRun]  # Total section width (mm)
        # Total section height (mm)
        self.t_height = data[0].loc['t_height', iRun]
        self.f_width = data[0].loc['f_width', iRun]  # Flange width (mm)
        self.f_thick = data[0].loc['f_thick', iRun]  # Flange thickness (mm)
        # Number of steel layers in flange ('1' or '2')
        self.f_lay = data[0].loc['f_lay', iRun]
        self.f_bar = data[0].loc['f_bar', iRun]  # Flange bar size (mm)
        self.f_spac = data[0].loc['f_spac', iRun]  # Flange bar spacing (mm)
        self.t_nt = data[0].loc['t_nt', iRun]  # Number of top bars
        self.t_tbar = data[0].loc['t_tbar', iRun]  # Diameter of top bars (mm)
        self.t_nb = data[0].loc['t_nb', iRun]  # Number of top bars
        self.t_bbar = data[0].loc['t_bbar', iRun]  # Diameter of top bars (mm)

        # Properties for Cruciform sections
        self.bX = data[0].loc['bX', iRun]  # Total section width (mm)
        self.bY = data[0].loc['bY', iRun]  # Total section height (mm)
        self.tX = data[0].loc['tX', iRun]  # Vertical flange thickness (mm)
        # Horizontal flange thickness (mm)
        self.tY = data[0].loc['tY', iRun]
        # Eccentricity of vertical flange (mm)
        self.cX = data[0].loc['cX', iRun]
        # Eccentricity of horizontal flange (mm)
        self.cY = data[0].loc['cY', iRun]
        # Diameter of shared bars in the center of the section (mm)
        self.cBar_centerD = data[0].loc['cBar_centerD', iRun]
        # Diameter of bars in the vertical flange (mm)
        self.cBar_vFlangeD = data[0].loc['cBar_vFlangeD', iRun]
        # Diameter of bars in the horizontal flange (mm)
        self.cBar_hFlangeD = data[0].loc['cBar_hFlangeD', iRun]
        # Number of bars in the vertical flange
        self.cBar_vFlangeN = data[0].loc['cBar_vFlangeN', iRun]
        # Number of bars in the horizontal flange
        self.cBar_hFlangeN = data[0].loc['cBar_hFlangeN', iRun]
        # Number of rows in the horizontal and vertical flanges
        self.cBar_rows = data[0].loc['cBar_rows', iRun]
        if self.section_type == 'Cruciform':
            self.cBar_rows.split(",")
            self.cBar_rows[0] = self.cBar_rows[0].replace('"', '')
            self.cBar_rows[1] = self.cBar_rows[1].replace('"', '')

        # Loading on section (M-φ diagram and shear calculations)
        # Applied moment in the section (kN-m)
        self.M_load = data[0].loc['M_load', iRun]
        # Applied axial force in the section (kN)
        self.N_load = data[0].loc['N_load', iRun]
        # Applied shear force in the section (kN)
        self.V_load = data[0].loc['V_load', iRun]

        # Shear reinforcement properties
        # Diameter of shear reinforcement (mm)
        self.d_shear = data[0].loc['d_shear', iRun]
        # Spacing of shear reinforcement (mm)
        self.s_shear = data[0].loc['s_shear', iRun]
        # Number of legs in section shear reinforcement
        self.n_legs = data[0].loc['n_legs', iRun].split(",")
        self.n_legs[0] = self.n_legs[0].replace('"', '')
        self.n_legs[1] = self.n_legs[1].replace('"', '')
        # Partial material safety factor for concrete
        self.gc_mat = data[0].loc['gc_mat', iRun]
        # Partial material safety factor for steel
        self.gs_mat = data[0].loc['gs_mat', iRun]
        # Design yield strength of shear reinforcement (N/mm2)
        self.Fywd = data[0].loc['Fywd', iRun]

        # Plastic hinge length
        # Plastic hinge length of the section in xx (mm or a factor)
        self.phl_xx = data[0].loc['phl_xx', iRun]
        # Plastic hinge length of the section in yy (mm or a factor)
        self.phl_yy = data[0].loc['phl_yy', iRun]

        # Parameters for failure criteria of rectangular beams
        # Dictated by property controlling failure (see ASCE 41-17 Table 10-17)
        self.rb_condition = data[0].loc['rb_condition', iRun]

        # Parameters for failure criteria of columns
        self.c_condition = data[0].loc['c_condition', iRun]
        self.ss_ratio_s = data[0].loc['ss_ratio_s', iRun]
        self.ss_ratio_t = data[0].loc['ss_ratio_t', iRun]

        # Handle area updating for initialization
        self.section_area = 0.0
        if self.section_type == 'Circular':
            self.section_area = 3.14 * \
                pow(float(self.diameter_circ) / 2000, 2.0)
        elif self.section_type == 'Rectangular Beam':
            self.section_area = float(
                self.rect_width) * float(self.rect_height) / (1000**2)
        elif self.section_type == 'Rectangular Column':
            self.section_area = round((float(self.rect_width) / 1000)
                                      * float(self.rect_height) / 1000, 5)
        elif self.section_type == 'T Beam':
            self.section_area = (float(self.f_width) * float(self.f_thick) +
                                 float(self.t_width) * (float(self.t_height) -
                                                        float(self.f_thick))) / (1000**2)
        elif self.section_type == 'Cruciform':
            # Check inputs
            if float(self.bX) - float(self.tX) <= 0 or float(self.bY) - \
                    float(self.tY) <= 0:
                raise AssertionError('Cruciform flange thickness is greater\
                                     than total section width!')
            if float(self.bX) - float(self.cX) <= 0 or float(self.bY) - \
                    float(self.cY) <= 0:
                raise AssertionError('Cruciform flange eccentricity is too\
                                     large!')
            self.section_area = round(((float(self.tX) *
                                        float(self.bY) +
                                        float(self.tY) *
                                        float(self.bX) -
                                        float(self.tX) *
                                        float(self.tY)) /
                                       (1000**2)), 5)

    def updateArea(self, area):
        # Handle area updating for logic loops
        if self.section_type == 'Circular':
            self.section_area = 3.14 * \
                pow(float(self.diameter_circ) / 2000, 2.0)
        elif self.section_type == 'Rectangular Beam':
            self.section_area = float(
                self.rect_width) * float(self.rect_height) / (1000**2)
        elif self.section_type == 'Rectangular Column':
            self.section_area = round((float(self.rect_width) / 1000)
                                      * float(self.rect_height) / 1000, 5)
        elif self.section_type == 'T Beam':
            self.section_area = (float(self.f_width) * float(self.f_thick) +
                                 float(self.t_width) * (float(self.t_height) -
                                                        float(self.f_thick))) / (1000**2)
        elif self.section_type == 'Cruciform':
            self.section_area = round(((float(self.tX) *
                                        float(self.bY) +
                                        float(self.tY) *
                                        float(self.bX) -
                                        float(self.tX) *
                                        float(self.tY)) /
                                       (1000**2)), 5)

        return self.section_area
