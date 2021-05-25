class info_elasticBeams():

    def __init__(self, data, iRun):

        ##################################################
        ### Static data (same for all batch instances) ###
        ##################################################

        # General project information
        self.author = data[1][0]  # Name of author
        self.company = data[1][1]  # Company name
        self.job_name = data[1][2]  # Project name

        # Batch mode options (False -> results only returned for given iSet)
        self.batchMode = data[1][3]
        # Input set for a single analysis (ignored if batchMode = True)
        self.iSet = data[1][4]
        self.iRun = int(iRun)
        self.tbOption = data[1][5]

        # Generate Unique Names for XTRACT
        # CAUTION: Modify standard values below with care
        # In almost all cases, editing is not required
        self.mat_name_unc = 'AutoUnconfined'  # Unconfinced concrete
        self.mat_name_conc = 'AutoConfined'  # Confined concrete
        self.mat_name_steel = 'AutoSteel'  # Steel
        self.section_name = 'AutoSection'  # Section name
        self.material_core = 'AutoConfined'  # Concrete core material name
        self.material_cover = 'AutoUnconfined'  # Concrete cover material name
        self.material_steel = 'AutoSteel'  # Steel material name

        # Specify path to save
        self.pSave = data[1][6]

        # Specify option to write terminal to log gile
        self.write_log = data[1][7]

        #####################
        ### Instance Data ###
        #####################

        # File name
        self.file_name = data[0].loc['info_filename', iRun]

        # Material data
        self.Ec = data[0].loc['Ec', iRun]  # Modulus of elasticity (MPa)
        self.den = data[0].loc['den', iRun]  # Density of concrete (kg/m3)
        self.v = data[0].loc['v', iRun]  # Poisson's ratio of concrete

        # General section properties
        # 'Rectangular Beam', 'Rectangular Column', 'Circular', 'T Beam'
        self.section_type = data[0].loc['section_type', iRun]  # Section type
        self.cover = data[0].loc['cover', iRun]  # concrete cover (mm)

        # Geometry for circular sections
        # Section Diameter (mm)
        self.diameter_circ = data[0].loc['diameter_circ', iRun]

        # Geometry for all rectangular sections
        self.rect_width = data[0].loc['rect_width', iRun]  # Section width (mm)
        # Section height (mm)
        self.rect_height = data[0].loc['rect_height', iRun]

        # Geometry for T Beam sections
        self.t_width = data[0].loc['t_width', iRun]  # Total section width (mm)
        # Total section height (mm)
        self.t_height = data[0].loc['t_height', iRun]
        self.f_width = data[0].loc['f_width', iRun]  # Flange width (mm)
        self.f_thick = data[0].loc['f_thick', iRun]  # Flange thickness (mm)

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

        # Section, material and curve nomenclature
        try:
            # Name sections and curves based on Beam ID
            self.beam_num = int(self.file_name[2:])

        except (IndexError, ValueError) as error:
            # Name sections and curves based on run number
            self.beam_num = int(self.iRun + 1) * 10

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
