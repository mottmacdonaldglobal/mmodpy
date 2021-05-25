# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 12:58:45 2018

@author: kevin.stanton
"""

import math


class Shear_Calcs():

    def __init__(self, aData):
        '''
        This constructor defines all the input parameters required to assess
        the shear capacity of a reinforced concrete section
        '''

        # Parameters for faliure criteria
        self.fParams = aData.fParams
        self.ss_ratio_s = float(aData.ss_ratio_s)
        self.ss_ratio_t = float(aData.ss_ratio_t)
        self.c_condition = aData.c_condition
        self.rb_condition = aData.rb_condition

        # Global properties
        self.c = float(aData.cover)  # Cover of reinforcement bars (mm)
        self.codeOption = aData.codeOption
        self.concreteFactor = aData.concreteFactor
        self.gs_mat = aData.gs_mat  # Partial factor for steel

        # Define section properties
        self.shape = aData.section_type
        if self.shape == 'Circular':

            # Diameter of circular section (mm)
            self.D = float(aData.diameter_circ)

            # Diameter of reinforcement bars (mm)
            self.db = float(aData.bar_diam_circ)

            # Number of reinforcement bars
            self.n = float(aData.num_bars_circ)

            # Diameter oustide stirrups (mm)
            self.Ds = self.D - 2 * self.c

        elif self.shape == 'Cruciform':

            # Equivalent section geometry for rectangular section components
            self.bx = float(aData.bX)  # loading in x
            self.by = float(aData.tY)  # loading in x
            self.hx = float(aData.bY)  # loading in y
            self.hy = float(aData.tX)  # loading in y

            # Average diameter of reinforcement bars in tension side
            self.db_tx = (float(aData.cBar_hFlangeD) +
                          float(aData.cBar_centerD)) / 2
            self.db_ty = (float(aData.cBar_vFlangeD) +
                          float(aData.cBar_centerD)) / 2

            # Average diameter of reinforcement bars in compression side
            self.db_cx = (float(aData.cBar_hFlangeD) +
                          float(aData.cBar_centerD)) / 2
            self.db_cy = (float(aData.cBar_vFlangeD) +
                          float(aData.cBar_centerD)) / 2

            # Number of reinforcement bars in tension side
            nBar_center = int(aData.cBar_rows[0]) * int(aData.cBar_rows[1])
            self.n_tx = int(aData.cBar_hFlangeN) / 2 + nBar_center / 2
            self.n_ty = int(aData.cBar_vFlangeN) / 2 + nBar_center / 2

            # Number of reinforcement bars in compression side
            self.n_cx = int(aData.cBar_hFlangeN) / 2 + nBar_center / 2
            self.n_cy = int(aData.cBar_vFlangeN) / 2 + nBar_center / 2

        elif (self.shape == 'Rectangular Beam') |\
             (self.shape == 'Rectangular Column'):

            # Dimensions of rectangular section for shear in x (mm)
            self.bx = float(aData.rect_height)  # width
            self.by = float(aData.rect_width)  # height

            # Dimensions of rectangular section for shear in y (mm)
            self.hx = float(aData.rect_width)  # width
            self.hy = float(aData.rect_height)  # height

            # Diameter of reinforcement bars in tension side (mm)
            if self.shape == 'Rectangular Beam':
                self.db_tx = float(aData.db_b)
                self.db_ty = (float(aData.db_b) + float(aData.db_t)) / 2

            elif self.shape == 'Rectangular Column':
                self.db_tx = float(aData.bar_diam_col)
                self.db_ty = float(aData.bar_diam_col)

            # Number of reinforcement bars in tension side
            if self.shape == 'Rectangular Beam':
                self.n_tx = float(aData.n_b)
                self.n_ty = float(2)  # Assumes no side bars
            elif self.shape == 'Rectangular Column':
                if int(aData.num_bars_col) == 8:
                    self.n_tx = 3
                    self.n_ty = 3
                elif int(aData.num_bars_col) == 4:
                    self.n_tx = 2
                    self.n_ty = 2
                elif int(aData.num_bars_col) == 6:
                    self.n_tx = 2
                    self.n_ty = 3
                elif int(aData.num_bars_col) == 10:
                    self.n_tx = 4
                    self.n_ty = 4
                elif int(aData.num_bars_col) == 12:
                    self.n_tx = 4
                    self.n_ty = 4
                elif int(aData.num_bars_col) == 14:
                    self.n_tx = 5
                    self.n_ty = 5

            # Diameter of reinforcement bars in compression side (mm)
            if self.shape == 'Rectangular Beam':
                self.db_cx = float(aData.db_t)
                self.db_cy = (float(aData.db_b) + float(aData.db_t)) / 2
            elif self.shape == 'Rectangular Column':
                self.db_cx = float(aData.bar_diam_col)
                self.db_cy = float(aData.bar_diam_col)

            # Number of reinforcement bars in compression side
            if self.shape == 'Rectangular Beam':
                self.n_cx = float(aData.n_t)
                self.n_cy = float(2)  # Assuming there are no side bars
            elif self.shape == 'Rectangular Column':
                if int(aData.num_bars_col) == 8:
                    self.n_cx = 3
                    self.n_cy = 3
                elif int(aData.num_bars_col) == 4:
                    self.n_cx = 2
                    self.n_cy = 2
                elif int(aData.num_bars_col) == 6:
                    self.n_cx = 2
                    self.n_cy = 3
                elif int(aData.num_bars_col) == 10:
                    self.n_cx = 4
                    self.n_cy = 4
                elif int(aData.num_bars_col) == 12:
                    self.n_cx = 4
                    self.n_cy = 4
                elif int(aData.num_bars_col) == 14:
                    self.n_cx = 5
                    self.n_cy = 5

        elif self.shape == 'T Beam':

            # Dimensions of T Beam section for loading in x (mm)
            self.bx = float(aData.t_height)  # width
            self.by = float(aData.t_width)  # height

            # Dimensions of T Beam section for loading in y (mm)
            self.hx = float(aData.t_width)  # width
            self.hy = float(aData.t_height)  # height

            # Diameter of reinforcement bars in tension side (mm)
            self.db_tx = float(aData.t_bbar)
            self.db_ty = (float(aData.t_bbar) + float(aData.t_tbar)) / 2

            # Number of reinforcement bars in tension side
            self.n_tx = float(aData.t_nb)
            self.n_ty = float(2)  # Assuming there are no side bars

            # Diameter of reinforcement bars in compression side (mm)
            self.db_cx = float(aData.t_tbar)
            self.db_cy = (float(aData.t_bbar) + float(aData.t_tbar)) / 2

            # Number of reinforcement bars in compression side
            self.n_cx = float(aData.t_nt)
            self.n_cy = float(2)  # Assuming there are no side bars

        # Diameter of shear reinforcement (mm)
        self.ds = float(aData.d_shear)

        # Spacing of shear reinforcement (mm)
        self.s = float(aData.s_shear)

        # Number of legs pertaining to shear reinforcement in the section
        self.nsx = int(aData.n_legs[0])
        self.nsy = int(aData.n_legs[1])

        # Define material properties
        # Modulus of elasticity of steel reinforcement (GPa)
        self.Est = float(aData.Es) / 1000

        # Characteristic yield strength of reinforcement (MPa)
        self.fyk = float(aData.Fy)

        # Partial material safety factor for concrete
        self.gc = float(aData.gc_mat)

        # Partial material safety factor for steel
        self.gs = float(aData.gs_mat)

        # Cylindrical characteristic compressive strength of concrete
        self.fck = float(aData.Fc)

        # Cylindrical characteristic mean strength of concrete
        self.fcm = float(aData.Fc)

        # Design yield strength of shear reinforcement (N/mm2)
        self.fywd = float(aData.Fywd)

        # Maximum unusable compressive strain of concrete
        self.ecu = float(aData.ecu)

        # Coefficient for long term effects on compressive strength
        self.acc = 1.0

        # Yield strain of steel reinforcement
        self.ey = self.fyk / self.Est

        # Design value of concrete compressive strength (MPa)
        self.fcd = self.acc * self.fcm / self.gc

        # Define loading conditions
        # Applied moment in the section (kN-m)
        self.Med = float(aData.M_load)

        # Applied axial force in the section (kN)
        self.Ned = float(aData.N_load)

        # Applied shear force in the section (kN)
        self.Ved = float(aData.V_load)
        self.file_name = aData.file_name

        # Define reinforcement details
        if self.shape == 'Circular':

            # Reinforcement area in tension side (mm2)
            self.As_t = self.n / 2 * math.pi * self.db**2 / 4

            # Reinforcement area in compression side (mm2)
            self.As_c = self.As_t

            # Area of nonprestressed reinforcement (mm2)
            self.Asc = self.As_t + self.As_c
            self.Ascx = self.Asc
            self.Ascy = self.Asc

            # Net area of section (mm2)
            self.Anc = math.pi * self.D**2 / 4 - self.Asc
            self.Ancx = self.Anc
            self.Ancy = self.Anc

            # Gross section area (mm2)
            self.Ag = math.pi * self.D**2 / 4
            self.Agcx = self.Ag
            self.Agcy = self.Ag

            # Effective depth in tension side (mm)
            self.d = 0.5 * self.D + 0.5 * self.Ds * math.tan(2 / math.pi)

            # Gross section area (mm2)
            self.AgShearcx = self.D * self.s
            self.AgShearcy = self.D * self.s

        elif (self.shape == 'Rectangular Beam') |\
             (self.shape == 'Rectangular Column') |\
             (self.shape == 'T Beam') |\
             (self.shape == 'Cruciform'):

            # Reinforcement area in tension side (mm2)
            self.As_tx = self.n_tx * math.pi * self.db_tx**2 / 4
            self.As_ty = self.n_ty * math.pi * self.db_ty**2 / 4

            # Reinforcement area in compression side (mm2)
            self.As_cx = self.n_cx * math.pi * self.db_cx**2 / 4
            self.As_cy = self.n_cy * math.pi * self.db_cy**2 / 4

            # Area of nonprestressed reinforcement (mm2)
            self.Ascx = self.As_tx + self.As_cx
            self.Ascy = self.As_ty + self.As_cy

            # Net area of section (mm2)
            self.Ancx = self.bx * self.by - self.Ascx
            self.Ancy = self.hx * self.hy - self.Ascy

            # Gross section area (mm2)
            self.Agcx = self.bx * self.by
            self.Agcy = self.hx * self.hy

            # Effective depth in tension side (mm)
            self.d_tx = self.bx - self.c - self.db_tx / 2
            self.d_ty = self.by - self.c - self.db_ty / 2

            # Effective depth in compression side (mm)
            self.d_cx = self.c + self.db_cx / 2
            self.d_cy = self.c + self.db_cy / 2

            # Gross section area (mm2)
            self.AgShearcx = self.bx * self.s
            self.AgShearcy = self.hx * self.s

    def __Conc_Cap__(self, Ned, width, height, n_t, d_t, n_b, d_b, As_t, As_b, As_total, A_conc_net, d_eff_tens):
        '''
        This function calculates the shear resistance provided by the concrete
        '''
        
        if self.codeOption == 'EC2':

            # Non-dimensional coefficient (EC2 §6.2.2)
            CRd_c = 0.18 / self.gc
            if (self.shape == 'Rectangular Beam') |\
                (self.shape == 'Rectangular Column') |\
                (self.shape == 'T Beam') |\
                    (self.shape == 'Cruciform'):

                # Non-dimensional coefficient (EC2 §6.2.2)
                k = min(2, 1 + math.sqrt(200 / d_eff_tens))

                # Non-dimensional coefficient (EC2 §6.2.2)
                k1 = 0.15

                # Smallest width of cross-section in the tensile area (mm)
                bw = min(width, height)

                # Minimum value of concrete shear strength (EC2 §6.2.2) (MPa)
                vmin = 0.035 * k**(3 / 2) * self.fcm**0.5

                # Ratio of tensile reinforcement (Asl)
                ρl = min(0.02, As_t / (bw * d_eff_tens))

                # Normal stress acting on the concrete section (MPa)
                σcp = min(0.2 * self.fcd, Ned / 1000 / (A_conc_net / 1000000))

                # Design shear resistance of the concrete alone (kN)
                VRd_c = max((CRd_c * k * (100 * ρl * self.fcm)**(1 / 3) +
                              k1 * σcp) * (bw * d_eff_tens) / 1000,
                             (vmin + k1 * σcp) * (bw * d_eff_tens) / 1000)

            elif self.shape == 'Circular':

                α = (2 / math.pi) * (self.Ds / self.D)

                # Non-dimensional coefficient (EC2 §6.2.2)
                Av = 1 / 4 * self.D**2 * \
                    (math.pi / 2 + α + math.sin(α) * math.cos(α))
                k = min(2, 1 + math.sqrt(200 / self.d))

                # Ratio of tensile reinforcement (Asl)
                ρl = As_t / (0.25 * math.pi * 1000**2)

                # Design shear resistance of the concrete alone (kN)
                VRd_c = (CRd_c * k * (100 * ρl * self.fcm)
                          ** (1 / 3)) * Av / 1000

        elif self.codeOption == 'CSA':

            ### CSA A23.3-14 ###

            # Coefficients
            beta = 0.18  # CONSTANT
            phiC = 1 / self.gc
            if self.concreteFactor == 'normal':
                cFac = 1.0
            elif self.concreteFactor == 'semi-light':
                cFac = 0.85
            elif self.concreteFactor == 'light':
                cFac = 0.75

            if (self.shape == 'Rectangular Beam') |\
                (self.shape == 'Rectangular Column') |\
                (self.shape == 'T Beam') |\
                    (self.shape == 'Cruciform'):

                # Effective shear depth
                dv = max(0.9 * d_eff_tens, 0.72 * height)

                # Smallest width of cross-section in the tensile area (mm)
                bw = width

                ε = min(max(1 /
                              (2 *
                               self.Est *
                               As_t) *
                              (abs(0.0) /
                               (dv /
                                1000) +
                                  0.0 -
                                  Ned *
                                  0.5), 0), 0.003)
                beta = 0.18
                phiC = 1 / self.gc

                # Design shear resistance of the concrete alone (kN)
                VRd_c = (phiC * cFac * beta * (self.fcm)
                          ** (0.5)) * (bw * dv) / 1E3

            elif self.shape == 'Circular':
                # Effective shear depth
                dv = max(0.9 * self.d, 0.72 * self.D)

                # Axial strain at Ned
                ε = min(max(1 /
                              (2 *
                               self.Est *
                               As_t) *
                              (abs(0.0) /
                               (dv /
                                1000) +
                                  0.0 -
                                  Ned *
                                  0.5), 0), 0.003)

                beta = 0.18
                # Design shear resistance of the concrete alone (kN)
                VRd_c = (phiC * cFac * beta * (self.fcm)
                          ** (0.5)) * self.d * self.D / 1E3

        # Return output (kN)
        return VRd_c

    def __Rein_Cap__(self, Ned, width, height, n_s, d_s, s_spacing, n_t, d_t, n_b, d_b, As_t, As_b, As_total, A_conc_net, d_eff_tens):
        '''
        Calculates the shear resistance provided by the shear
        reinforcement
        '''
        
        # Cross-sectional area of the shear reinforcement (mm2)
        Asw = n_s * math.pi * (d_s**2) / 4
        self.Aswx = Asw # for legacy fParams logic
        self.Aswy = Asw # for legacy fParams logic

        if self.codeOption == 'EC2':

            # Strength reduction factor for concrete cracked in shear (EC2 §6.2.2)
            v = 0.6 * (1 - self.fcm / 250)

            if (self.shape == 'Rectangular Beam') |\
                (self.shape == 'Rectangular Column') |\
                (self.shape == 'T Beam') |\
                    (self.shape == 'Cruciform'):

                # Inner lever arm, for a member with constant depth (mm)
                z = 0.9 * d_eff_tens

                # Longitudinal strain at mid-depth (MC2010, 7.3-16)
                ε = min(max((1 / (2 * self.Est * As_t)) *
                              (abs(self.Med) / (z / 1000) + self.Ved - Ned * 0.5), 0), 0.003)

                # Αngle of concrete compression strut (MC2010, 7.3-39)
                # (degrees)
                θmin = 20 + 10000 * ε

                # Angle of Shear Reinforcement to Beam Axis (degrees)
                alpha = 90

                # Coefficient considering the state of stress in compression
                # chord
                acw = 1

                # Smallest width of cross-section in the tensile area (mm)
                bw = min(width, height)

                # Shear resistance due to crushing of compression struts (kN)
                VRd_max = abs(acw *
                               bw *
                               z *
                               v *
                               self.fcd /
                               (1000 *
                                (1 /
                                 math.tan(math. radians(θmin)) +
                                    1 /
                                    math.tan(math.radians(alpha)))))

                # Design shear resistance provided by shear reinforcement (kN)
                VRd_s = min(abs(Asw /
                                 s_pacing *
                                 z *
                                 self.fywd *
                                 (1 /
                                  math.tan(math.radians(θmin)) +
                                     1 /
                                     math.tan(math.radians(alpha))) /
                                 1000), VRd_max)

            elif self.shape == 'Circular':

                # Compute alpha
                α = (2 / math.pi) * (self.Ds / self.D)

                # Inner lever arm, for a member with constant depth (mm)
                z = max(0.9 * self.d, 0.72 * self.D)

                # Longitudinal strain at mid-depth (MC2010, 7.3-16)
                ε = min(max(1 / (2 * self.Est * As_t) * (self.Med / \
                          (z / 1000) + self.Ved - Ned * 0.5), 0), 0.003) # m/m

                # Αngle of concrete compression strut (MC2010, 7.3-39)
                # (degrees)
                θmin = 20 + 10000 * ε

                # Non-dimensional coefficient (EC2 §6.2.2)
                Av = 1 / 4 * self.D**2 * \
                    (math.pi / 2 + α + math.sin(α) * math.cos(α))

                # Shear resistance due to crushing of compression struts (kN)
                VRd_max = abs(v *
                               self.fcd *
                               Av *
                               z /
                               (1 /
                                math.tan(math.radians(θmin)) +
                                   math.tan(math.radians(θmin))) /
                               1000)                

                # Design shear resistance provided by shear reinforcement (kN)
                VRd_s = min(1.34 * Asw * self.fywd / s_spacing * self.D *
                             (1 / math.tan(math. radians(θmin))) / 1000, VRd_max)

        if self.codeOption == 'CSA':

            phiC = 1 / self.gc

            if (self.shape == 'Rectangular Beam') |\
                (self.shape == 'Rectangular Column') |\
                (self.shape == 'T Beam') |\
                    (self.shape == 'Cruciform'):

                # Inner lever arm, for a member with constant depth (mm)
                dv = max(0.9 * d_eff_tens, 0.72 * height)
                # Longitudinal strain at mid-depth
                ε = min(max(1 /
                              (2 *
                               self.Est *
                               As_t) *
                              (abs(0.0) /
                               (dv /
                                1000) +
                                  0.0 -
                                  Ned *
                                  0.5), 0), 0.003)

                # Αngle of concrete compression strut (degrees)
                θmin = 29 + 7000 * ε

                # Smallest width of cross-section in the tensile area (mm)
                bw = min(width, height)

                # Design shear resistance provided by shear reinforcement (kN)
                VRd_s = phiC * Asw * self.fywd * dv * \
                    (1 / math.tan(math.radians(θmin))) / (1000 * s_spacing)

            elif self.shape == 'Circular':

                # Compute alpha
                α = (2 / math.pi) * (self.Ds / self.D)

                # Inner lever arm, for a member with constant depth (mm)
                z = max(0.9 * self.d, 0.72 * self.D)

                # Longitudinal strain at mid-depth (MC2010, 7.3-16)
                ε = min(max(1 / (2 * self.Est * As_t) *
                              (0.0 / (z / 1000) + 0.0 - Ned * 0.5), 0), 0.003)

                # Αngle of concrete compression strut (MC2010, 7.3-39)
                # (degrees)
                θmin = 29 + 7000 * ε

                # Non-dimensional coefficient (EC2 §6.2.2)
                Av = 1 / 4 * self.D**2 * \
                    (math.pi / 2 + α + math.sin(α) * math.cos(α))

                # Design shear resistance provided by shear reinforcement (kN)
                VRd_s = phiC * Asw * self.fywd * z * \
                    (1 / math.tan(math.radians(θmin))) / (1000 * s_spacing)

        # Return output (kN)
        return VRd_s

    def Shear_cap(self, Ned, width, height, n_s, d_s, s_spacing, n_t, d_t, n_b, d_b):
        '''
        Calculates shear capacity of a reinforced concrete section
        '''
        
        # Area of steel in top/bottom (tension/compression)
        As_t = n_t * math.pi * (d_t/2)**2
        As_b = n_b * math.pi * (d_b/2)**2
        # Total steel area (mm2)
        As_total = As_t + As_b
        # Net concrete area of section (mm2)
        A_conc_net = width * height - As_total
        # Effective depth in tension side (mm)
        d_eff_tens = height - self.c - d_t / 2

        # Calculate shear capacity of concrete
        VRd_c = self.__Conc_Cap__(Ned, width, height, n_t, d_t, n_b, d_b, As_t, As_b, As_total, A_conc_net, d_eff_tens)

        # Calculate shear capacity of reinforcement
        VRd_s = self.__Rein_Cap__(Ned, width, height, n_s, d_s, s_spacing, n_t, d_t, n_b, d_b, As_t, As_b, As_total, A_conc_net, d_eff_tens)

        # Define shear capacity of reinforced concrete section
        if self.codeOption == 'EC2':
            VRd = max(VRd_c, VRd_s)
        if self.codeOption == 'CSA':
            phiC = 1 / self.gc
            if self.shape == 'Circular':
                dv = self.d
                bw = self.D
            else:
                dv = max(0.9 * d_eff_tens, 0.72 * height)
                bw = width
            VRd = min(VRd_c + VRd_s, 0.25 * phiC *
                       self.fcm * (bw * dv) / 1E3)

        # Return output
        return VRd

    def interpLin(self, y0, y1, x0, x1, x):
        y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
        return y

    def prLimits(self):
        '''
        Calculates plastic rotation thresholds following procedures from ASCE 41
        '''

        if self.fParams:
            if (self.Ved == 0) |\
                (self.Ved == 0) |\
                (math.isnan(self.Ved)):
                print(
                    'Warning: zero shear demand assumed (check inputs if unexpected)')
                self.Ved = 0.001
                
            if (self.shape == 'Rectangular Beam') |\
                    (self.shape == 'T Beam'):   
                # ASCE 41-17 Table 10-7
                # Note: ASCE 41-17 is the same as ASCE 41-13 for beams                
                if any("1" in i for i in self.rb_condition):

                    # Ratio of nonprestressed tension/compression reinforcement
                    # to net area for s/t-directions
                    rhoTs = self.As_tx / self.Ancx
                    rhoTt = self.As_ty / self.Ancy
                    rhoCs = self.As_cx / self.Ancx
                    rhoCt = self.As_cy / self.Ancy

                    # Compute reinforcement ratio producing balanced strain
                    if self.fck <= 30:
                        beta = 0.85
                    else:
                        beta = max(0.85 - 0.05 / (7 * (self.fck - 30)), 0.65)
                    self.rho_bal_s = 0.85 * self.fck * beta * self.ecu * self.Est * \
                        0.001 / (self.fyk * (self.ecu * self.Est * 0.001 + self.fyk))
                    self.rho_bal_t = self.rho_bal_s

                    # Determine balance point in s/t-directions
                    bP_s = (rhoTs - rhoCs) / self.rho_bal_s
                    bP_t = (rhoTt - rhoCt) / self.rho_bal_t

                    # Determine compliance of transverse reinforcement
                    self.b = max(self.d_tx, self.d_ty)  # Governing dimension
                    if self.s <= self.b / 3:
                        compliance_t = 'C'
                        compliance_s = 'C'
                    else:
                        compliance_s = 'NC'
                        compliance_t = 'NC'

                    # Determine shear constant
                    vC_s = self.Ved * 1000 / \
                        (self.by * self.bx * 0.8 * (self.fck)**0.5)
                    vC_t = self.Ved * 1000 / \
                        (self.bx * self.by * 0.8 * (self.fck)**0.5)

                    # Determine rotation limits defining failure in the local
                    # s-direction based on Table 10-7 (ASCE 41-17)
                    if compliance_s == 'C':
                        if vC_s <= 0.25 and bP_s <= 0:
                            # If row 1
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.01, 0.025, 0.05
                        elif vC_s >= 0.25 and vC_s <= 0.5 and bP_s <= 0:
                            # If between row 1 and 2 for vC_s
                            c1_PRS1 = 0
                            c1_PRS2 = self.interpLin(
                                0.01, 0.005, 0.25, 0.5, vC_s)
                            c1_PRS3 = self.interpLin(
                                0.025, 0.02, 0.25, 0.5, vC_s)
                            c1_PRS4 = self.interpLin(
                                0.05, 0.04, 0.25, 0.5, vC_s)
                        elif vC_s >= 0.5 and bP_s <= 0:
                            # If row 2
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.005, 0.02, 0.04
                        elif vC_s <= 0.25 and bP_s >= 0 and bP_s <= 0.5:
                            # If between row 1 and 3 for bP_s
                            c1_PRS1 = 0
                            c1_PRS2 = self.interpLin(0.01, 0.005, 0, 0.5, bP_s)
                            c1_PRS3 = self.interpLin(0.025, 0.02, 0, 0.5, bP_s)
                            c1_PRS4 = self.interpLin(0.05, 0.03, 0, 0.5, bP_s)
                        elif vC_s <= 0.25 and bP_s >= 0.5:
                            # If row 3
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.005, 0.02, 0.03
                        elif vC_s >= 0.5 and bP_s >= 0 and bP_s <= 0.5:
                            # If between row 2 and 4 for bP_s
                            c1_PRS1 = 0
                            c1_PRS2 = 0.005
                            c1_PRS3 = self.interpLin(0.02, 0.015, 0, 0.5, bP_s)
                            c1_PRS4 = self.interpLin(0.04, 0.04, 0, 0.5, bP_s)
                        elif vC_s >= 0.25 and vC_s <= 0.5 and bP_s >= 0.5:
                            # If between row 3 and 4 for vC_s
                            c1_PRS1 = 0
                            c1_PRS2 = 0.005
                            c1_PRS3 = self.interpLin(
                                0.02, 0.015, 0.25, 0.5, bP_s)
                            c1_PRS4 = self.interpLin(
                                0.03, 0.02, 0.25, 0.5, bP_s)
                        elif vC_s >= 0.5 and bP_s >= 0.5:
                            # If row 4
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.005, 0.015, 0.02

                    elif compliance_s == 'NC':
                        if vC_s <= 0.25 and bP_s <= 0:
                            # If row 1
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.005, 0.02, 0.03
                        elif vC_s >= 0.25 and vC_s <= 0.5 and bP_s <= 0:
                            # If between row 1 and 2 for vC_s
                            c1_PRS1 = 0
                            c1_PRS2 = self.interpLin(
                                0.005, 0.0015, 0, 0.5, bP_s)
                            c1_PRS3 = self.interpLin(0.02, 0.01, 0, 0.5, bP_s)
                            c1_PRS4 = self.interpLin(0.03, 0.015, 0, 0.5, bP_s)
                        elif vC_s >= 0.5 and bP_s <= 0:
                            # If row 2
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.0015, 0.01, 0.015
                        elif vC_s <= 0.25 and bP_s >= 0 and bP_s <= 0.5:
                            # If between row 1 and 3 for bP_s
                            c1_PRS1 = 0
                            c1_PRS2 = 0.005
                            c1_PRS3 = self.interpLin(0.02, 0.01, 0, 0.5, bP_s)
                            c1_PRS4 = self.interpLin(0.03, 0.015, 0, 0.5, bP_s)
                        elif vC_s <= 0.25 and bP_s >= 0.5:
                            # If row 3
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.005, 0.01, 0.015
                        elif vC_s >= 0.5 and bP_s >= 0 and bP_s <= 0.5:
                            # If between row 2 and 4 for bP_s
                            c1_PRS1 = 0
                            c1_PRS2 = 0.0015
                            c1_PRS3 = self.interpLin(0.01, 0.005, 0, 0.5, bP_s)
                            c1_PRS4 = self.interpLin(0.015, 0.01, 0, 0.5, bP_s)
                        elif vC_s >= 0.25 and vC_s <= 0.5 and bP_s >= 0.5:
                            # If between row 3 and 4 for vC_s
                            c1_PRS1 = 0
                            c1_PRS2 = self.interpLin(
                                0.005, 0.0015, 0.25, 0.5, bP_s)
                            c1_PRS3 = self.interpLin(
                                0.01, 0.005, 0.25, 0.5, bP_s)
                            c1_PRS4 = self.interpLin(
                                0.015, 0.01, 0.25, 0.5, bP_s)
                        elif vC_s >= 0.5 and bP_s >= 0.5:
                            # If row 4
                            c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.0015, 0.005, 0.01

                    if compliance_t == 'C':
                        if vC_t <= 0.25 and bP_t <= 0:
                            # If row 1
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.01, 0.025, 0.05
                        elif vC_t >= 0.25 and vC_t <= 0.5 and bP_t <= 0:
                            # If between row 1 and 2 for vC_t
                            c1_PRT1 = 0
                            c1_PRT2 = self.interpLin(
                                0.01, 0.005, 0.25, 0.5, vC_t)
                            c1_PRT3 = self.interpLin(
                                0.025, 0.02, 0.25, 0.5, vC_t)
                            c1_PRT4 = self.interpLin(
                                0.05, 0.04, 0.25, 0.5, vC_t)
                        elif vC_t >= 0.5 and bP_t <= 0:
                            # If row 2
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.005, 0.02, 0.04
                        elif vC_t <= 0.25 and bP_t >= 0 and bP_t <= 0.5:
                            # If between row 1 and 3 for bP_t
                            c1_PRT1 = 0
                            c1_PRT2 = self.interpLin(0.01, 0.005, 0, 0.5, bP_t)
                            c1_PRT3 = self.interpLin(0.025, 0.02, 0, 0.5, bP_t)
                            c1_PRT4 = self.interpLin(0.05, 0.03, 0, 0.5, bP_t)
                        elif vC_t <= 0.25 and bP_t >= 0.5:
                            # If row 3
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.005, 0.02, 0.03
                        elif vC_t >= 0.5 and bP_t >= 0 and bP_t <= 0.5:
                            # If between row 2 and 4 for bP_t
                            c1_PRT1 = 0
                            c1_PRT2 = 0.005
                            c1_PRT3 = self.interpLin(0.02, 0.015, 0, 0.5, bP_t)
                            c1_PRT4 = self.interpLin(0.04, 0.04, 0, 0.5, bP_t)
                        elif vC_t >= 0.25 and vC_t <= 0.5 and bP_t >= 0.5:
                            # If between row 3 and 4 for vC_t
                            c1_PRT1 = 0
                            c1_PRT2 = 0.005
                            c1_PRT3 = self.interpLin(
                                0.02, 0.015, 0.25, 0.5, bP_t)
                            c1_PRT4 = self.interpLin(
                                0.03, 0.02, 0.25, 0.5, bP_t)
                        elif vC_t >= 0.5 and bP_t >= 0.5:
                            # If row 4
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.005, 0.015, 0.02

                    elif compliance_t == 'NC':
                        if vC_t <= 0.25 and bP_t <= 0:
                            # If row 1
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.005, 0.02, 0.03
                        elif vC_t >= 0.25 and vC_t <= 0.5 and bP_t <= 0:
                            # If between row 1 and 2 for vC_t
                            c1_PRT1 = 0
                            c1_PRT2 = self.interpLin(
                                0.005, 0.0015, 0, 0.5, bP_t)
                            c1_PRT3 = self.interpLin(0.02, 0.01, 0, 0.5, bP_t)
                            c1_PRT4 = self.interpLin(0.03, 0.015, 0, 0.5, bP_t)
                        elif vC_t >= 0.5 and bP_t <= 0:
                            # If row 2
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.0015, 0.01, 0.015
                        elif vC_t <= 0.25 and bP_t >= 0 and bP_t <= 0.5:
                            # If between row 1 and 3 for bP_t
                            c1_PRT1 = 0
                            c1_PRT2 = 0.005
                            c1_PRT3 = self.interpLin(0.02, 0.01, 0, 0.5, bP_t)
                            c1_PRT4 = self.interpLin(0.03, 0.015, 0, 0.5, bP_t)
                        elif vC_t <= 0.25 and bP_t >= 0.5:
                            # If row 3
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.005, 0.01, 0.015
                        elif vC_t >= 0.5 and bP_t >= 0 and bP_t <= 0.5:
                            # If between row 2 and 4 for bP_t
                            c1_PRT1 = 0
                            c1_PRT2 = 0.0015
                            c1_PRT3 = self.interpLin(0.01, 0.005, 0, 0.5, bP_t)
                            c1_PRT4 = self.interpLin(0.015, 0.01, 0, 0.5, bP_t)
                        elif vC_t >= 0.25 and vC_t <= 0.5 and bP_t >= 0.5:
                            # If between row 3 and 4 for vC_t
                            c1_PRT1 = 0
                            c1_PRT2 = self.interpLin(
                                0.005, 0.0015, 0.25, 0.5, bP_t)
                            c1_PRT3 = self.interpLin(
                                0.01, 0.005, 0.25, 0.5, bP_t)
                            c1_PRT4 = self.interpLin(
                                0.015, 0.01, 0.25, 0.5, bP_t)
                        elif vC_t >= 0.5 and bP_t >= 0.5:
                            # If row 4
                            c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.0015, 0.005, 0.01

                if any("2" in i for i in self.rb_condition):
                    if self.s <= self.bx / 2:
                        c2_PRS1, c2_PRS2, c2_PRS3, c2_PRS4 = 0, 0.0015, 0.01, 0.02
                    else:
                        c2_PRS1, c2_PRS2, c2_PRS3, c2_PRS4 = 0, 0.0015, 0.005, 0.01
                    if self.s <= self.by / 2:
                        c2_PRT1, c2_PRT2, c2_PRT3, c2_PRT4 = 0, 0.0015, 0.01, 0.02
                    else:
                        c2_PRT1, c2_PRT2, c2_PRT3, c2_PRT4 = 0, 0.0015, 0.005, 0.01

                if any("3" in i for i in self.rb_condition):
                    if self.s <= self.bx / 2:
                        c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 0, 0.0015, 0.01, 0.02
                    else:
                        c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 0, 0.0015, 0.005, 0.01
                    if self.s <= self.by / 2:
                        c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 0, 0.0015, 0.01, 0.02
                    else:
                        c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 0, 0.0015, 0.005, 0.01

                if any("4" in i for i in self.rb_condition):
                    c4_PRS1, c4_PRS2, c4_PRS3, c4_PRS4 = 0, 0.01, 0.02, 0.03
                    c4_PRT1, c4_PRT2, c4_PRT3, c4_PRT4 = 0, 0.01, 0.02, 0.03

                # Specify dummy values for conditions that do not control
                if not any("1" in i for i in self.rb_condition):
                    c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 1, 1, 1, 1
                    c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 1, 1, 1, 1
                if not any("2" in i for i in self.rb_condition):
                    c2_PRS1, c2_PRS2, c2_PRS3, c2_PRS4 = 1, 1, 1, 1
                    c2_PRT1, c2_PRT2, c2_PRT3, c2_PRT4 = 1, 1, 1, 1
                if not any("3" in i for i in self.rb_condition):
                    c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 1, 1, 1, 1
                    c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 1, 1, 1, 1
                if not any("4" in i for i in self.rb_condition):
                    c4_PRS1, c4_PRS2, c4_PRS3, c4_PRS4 = 1, 1, 1, 1
                    c4_PRT1, c4_PRT2, c4_PRT3, c4_PRT4 = 1, 1, 1, 1

                # Take the minimum outcomes among all specified conditions
                PRS1, PRS2, PRS3, PRS4 = [
                    min(
                        c1_PRS1, c2_PRS1, c3_PRS1, c4_PRS1), min(
                        c1_PRS2, c2_PRS2, c3_PRS2, c4_PRS2), min(
                        c1_PRS3, c2_PRS3, c3_PRS3, c4_PRS3), min(
                        c1_PRS4, c2_PRS4, c3_PRS4, c4_PRS4)]
                PRT1, PRT2, PRT3, PRT4 = [
                    min(
                        c1_PRT1, c2_PRT1, c3_PRT1, c4_PRT1), min(
                        c1_PRT2, c2_PRT2, c3_PRT2, c4_PRT2), min(
                        c1_PRT3, c2_PRT3, c3_PRT3, c4_PRT3), min(
                        c1_PRT4, c2_PRT4, c3_PRT4, c4_PRT4)]

            if (self.shape == 'Circular') |\
                    (self.shape == 'Rectangular Column') |\
                    (self.shape == 'Cruciform'):
                
                if self.fParams == 'ASCE 41-13':
                    # General parameters for ASCE 41-13 Table 10-8
                    paf = self.Ned/(self.Ag*self.fck)
                    rho_s = self.AgShearcx/(self.bx*self.s)
                    rho_t = self.AgShearcy/(self.by*self.s)
                    # Determine shear constant
                    vC_s = self.Ved * 1000 / \
                        (self.by * 0.8 * self.bx * (self.fck)**0.5)
                    vC_t = self.Ved * 1000 / \
                        (self.bx * 0.8 * self.by * (self.fck)**0.5)
                    
                    if any("1" in i for i in self.c_condition):
                        c1_PRS1 = 0
                        c1_PRT1 = 0
                        if paf <= 0.1:
                            if rho_s >= 0.006:
                                c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.005, 0.045, 0.06
                            elif rho_s <= 0.002:
                                c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.005, 0.027, 0.034
                            else:
                                c1_PRS2 = 0.005
                                c1_PRS3 = self.interpLin(0.027, 0.045, 0.002, 0.006, rho_s)
                                c1_PRS4 = self.interpLin(0.034, 0.06, 0.002, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.005, 0.045, 0.06                                
                            elif rho_t <= 0.002:
                                c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.005, 0.027, 0.034
                            else:
                                c1_PRT2 = 0.005
                                c1_PRT3 = self.interpLin(0.027, 0.045, 0.002, 0.006, rho_t)
                                c1_PRT4 = self.interpLin(0.034, 0.06, 0.002, 0.006, rho_t)
                        if paf >= 0.6:
                            if rho_s >= 0.006:
                                c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.003, 0.009, 0.01
                            elif rho_s <= 0.002:
                                c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 0, 0.002, 0.004, 0.005
                            else:
                                c1_PRS2 = self.interpLin(0.002, 0.003, 0.002, 0.006, rho_s)
                                c1_PRS3 = self.interpLin(0.004, 0.009, 0.002, 0.006, rho_s)
                                c1_PRS4 = self.interpLin(0.005, 0.01, 0.002, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.003, 0.009, 0.01                               
                            elif rho_t <= 0.002:
                                c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 0, 0.002, 0.004, 0.005
                            else:
                                c1_PRT2 = self.interpLin(0.002, 0.003, 0.002, 0.006, rho_t)
                                c1_PRT3 = self.interpLin(0.004, 0.009, 0.002, 0.006, rho_t)
                                c1_PRT4 = self.interpLin(0.005, 0.01, 0.002, 0.006, rho_t)
                        else:
                            c1_PRS2_g = self.interpLin(0.003, 0.005, 0.6, 0.1, paf)
                            c1_PRS3_g = self.interpLin(0.009, 0.045, 0.6, 0.1, paf)
                            c1_PRS4_g = self.interpLin(0.01, 0.06, 0.6, 0.1, paf)
                            c1_PRS2_l = self.interpLin(0.002, 0.005, 0.6, 0.1, paf)
                            c1_PRS3_l = self.interpLin(0.004, 0.027, 0.6, 0.1, paf)
                            c1_PRS4_l = self.interpLin(0.005, 0.034, 0.6, 0.1, paf)
                            c1_PRT2_g = self.interpLin(0.003, 0.005, 0.6, 0.1, paf)
                            c1_PRT3_g = self.interpLin(0.009, 0.045, 0.6, 0.1, paf)
                            c1_PRT4_g = self.interpLin(0.01, 0.06, 0.6, 0.1, paf)
                            c1_PRT2_l = self.interpLin(0.002, 0.005, 0.6, 0.1, paf)
                            c1_PRT3_l = self.interpLin(0.004, 0.027, 0.6, 0.1, paf)
                            c1_PRT4_l = self.interpLin(0.005, 0.034, 0.6, 0.1, paf)
                            if rho_s >= 0.006:
                                c1_PRS2 = c1_PRS2_g
                                c1_PRS3 = c1_PRS3_g
                                c1_PRS4 = c1_PRS4_g
                            elif rho_s <= 0.002:
                                c1_PRS2 = c1_PRS2_l
                                c1_PRS3 = c1_PRS3_l
                                c1_PRS4 = c1_PRS4_l
                            else:
                                c1_PRS2 = self.interpLin(c1_PRS2_l, c1_PRS2_g, 0.002, 0.006, rho_s)
                                c1_PRS3 = self.interpLin(c1_PRS3_l, c1_PRS3_g, 0.002, 0.006, rho_s)
                                c1_PRS4 = self.interpLin(c1_PRS4_l, c1_PRS4_g, 0.002, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c1_PRT2 = c1_PRT2_g
                                c1_PRT3 = c1_PRT3_g
                                c1_PRT4 = c1_PRT4_g
                            elif rho_t <= 0.002:
                                c1_PRT2 = c1_PRT2_l
                                c1_PRT3 = c1_PRT3_l
                                c1_PRT4 = c1_PRT4_l
                            else:
                                c1_PRT2 = self.interpLin(c1_PRT2_l, c1_PRT2_g, 0.002, 0.006, rho_t)
                                c1_PRT3 = self.interpLin(c1_PRT3_l, c1_PRT3_g, 0.002, 0.006, rho_t)
                                c1_PRT4 = self.interpLin(c1_PRT4_l, c1_PRT4_g, 0.002, 0.006, rho_t)
                                
                    if any("2" in i for i in self.c_condition):
                        c2_PRS1 = 0
                        c2_PRT1 = 0
                        if paf <= 0.1:
                            if rho_s >= 0.006:
                                c2_PRS1_vl, c2_PRS2_vl, c2_PRS3_vl, c2_PRS4_vl = 0, 0.005, 0.045, 0.06
                                c2_PRS1_vg, c2_PRS2_vg, c2_PRS3_vg, c2_PRS4_vg = 0, 0.005, 0.045, 0.06
                            elif rho_s <= 0.0005:
                                c2_PRS1_vl, c2_PRS2_vl, c2_PRS3_vl, c2_PRS4_vl = 0, 0.005, 0.01, 0.012
                                c2_PRS1_vg, c2_PRS2_vg, c2_PRS3_vg, c2_PRS4_vg = 0, 0.004, 0.005, 0.006
                            else:
                                c2_PRS2_vl = 0.005
                                c2_PRS3_vl = self.interpLin(0.01, 0.045, 0.0005, 0.006, rho_s)
                                c2_PRS4_vl = self.interpLin(0.012, 0.06, 0.0005, 0.006, rho_s)     
                                c2_PRS2_vg = self.interpLin(0.004, 0.005, 0.0005, 0.006, rho_s)
                                c2_PRS3_vg = self.interpLin(0.005, 0.045, 0.0005, 0.006, rho_s)
                                c2_PRS4_vg = self.interpLin(0.006, 0.06, 0.0005, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c2_PRT1_vl, c2_PRT2_vl, c2_PRT3_vl, c2_PRT4_vl = 0, 0.005, 0.045, 0.06
                                c2_PRT1_vg, c2_PRT2_vg, c2_PRT3_vg, c2_PRT4_vg = 0, 0.005, 0.045, 0.06
                            elif rho_t <= 0.0005:
                                c2_PRT1_vl, c2_PRT2_vl, c2_PRT3_vl, c2_PRT4_vl = 0, 0.005, 0.01, 0.012
                                c2_PRT1_vg, c2_PRT2_vg, c2_PRT3_vg, c2_PRT4_vg = 0, 0.004, 0.005, 0.006
                            else:
                                c2_PRT2_vl = 0.005
                                c2_PRT3_vl = self.interpLin(0.01, 0.045, 0.0005, 0.006, rho_t)
                                c2_PRT4_vl = self.interpLin(0.012, 0.06, 0.0005, 0.006, rho_t)     
                                c2_PRT2_vg = self.interpLin(0.004, 0.005, 0.0005, 0.006, rho_t)
                                c2_PRT3_vg = self.interpLin(0.005, 0.045, 0.0005, 0.006, rho_t)
                                c2_PRT4_vg = self.interpLin(0.006, 0.06, 0.0005, 0.006, rho_t)
                        if paf >= 0.6:
                            if rho_s >= 0.006:
                                c2_PRS1_vl, c2_PRS2_vl, c2_PRS3_vl, c2_PRS4_vl = 0, 0.003, 0.009, 0.01
                                c2_PRS1_vg, c2_PRS2_vg, c2_PRS3_vg, c2_PRS4_vg = 0, 0.003, 0.007, 0.008
                            elif rho_s <= 0.0005:
                                c2_PRS1_vl, c2_PRS2_vl, c2_PRS3_vl, c2_PRS4_vl = 0, 0.002, 0.003, 0.004
                                c2_PRS1_vg, c2_PRS2_vg, c2_PRS3_vg, c2_PRS4_vg = 0, 0.000001, 0.000001, 0.000001
                            else:
                                c2_PRS2_vl = self.interpLin(0.002, 0.003, 0.0005, 0.006, rho_s)
                                c2_PRS3_vl = self.interpLin(0.003, 0.009, 0.0005, 0.006, rho_s)
                                c2_PRS4_vl = self.interpLin(0.004, 0.01, 0.0005, 0.006, rho_s) 
                                c2_PRS2_vg = self.interpLin(0.000001, 0.003, 0.0005, 0.006, rho_s)
                                c2_PRS3_vg = self.interpLin(0.000001, 0.007, 0.0005, 0.006, rho_s)
                                c2_PRS4_vg = self.interpLin(0.000001, 0.008, 0.0005, 0.006, rho_s)   
                            if rho_t >= 0.006:
                                c2_PRT1_vl, c2_PRT2_vl, c2_PRT3_vl, c2_PRT4_vl = 0, 0.003, 0.009, 0.01
                                c2_PRT1_vg, c2_PRT2_vg, c2_PRT3_vg, c2_PRT4_vg = 0, 0.003, 0.007, 0.008
                            elif rho_t <= 0.0005:
                                c2_PRT1_vl, c2_PRT2_vl, c2_PRT3_vl, c2_PRT4_vl = 0, 0.002, 0.003, 0.004
                                c2_PRT1_vg, c2_PRT2_vg, c2_PRT3_vg, c2_PRT4_vg = 0, 0.000001, 0.000001, 0.000001
                            else:
                                c2_PRT2_vl = self.interpLin(0.002, 0.003, 0.0005, 0.006, rho_t)
                                c2_PRT3_vl = self.interpLin(0.003, 0.009, 0.0005, 0.006, rho_t)
                                c2_PRT4_vl = self.interpLin(0.004, 0.01, 0.0005, 0.006, rho_t) 
                                c2_PRT2_vg = self.interpLin(0.000001, 0.003, 0.0005, 0.006, rho_t)
                                c2_PRT3_vg = self.interpLin(0.000001, 0.007, 0.0005, 0.006, rho_t)
                                c2_PRT4_vg = self.interpLin(0.000001, 0.008, 0.0005, 0.006, rho_t)
                        else:
                            c2_PRS2_g_vg = self.interpLin(0.003, 0.005, 0.6, 0.1, paf)
                            c2_PRS3_g_vg = self.interpLin(0.007, 0.045, 0.6, 0.1, paf)
                            c2_PRS4_g_vg = self.interpLin(0.008, 0.06, 0.6, 0.1, paf)
                            c2_PRS2_g_vl = self.interpLin(0.003, 0.005, 0.6, 0.1, paf)
                            c2_PRS3_g_vl = self.interpLin(0.009, 0.045, 0.6, 0.1, paf)
                            c2_PRS4_g_vl = self.interpLin(0.01, 0.06, 0.6, 0.1, paf)
                            c2_PRS2_l_vg = self.interpLin(0.000001, 0.004, 0.6, 0.1, paf)
                            c2_PRS3_l_vg = self.interpLin(0.000001, 0.005, 0.6, 0.1, paf)
                            c2_PRS4_l_vg = self.interpLin(0.000001, 0.006, 0.6, 0.1, paf)
                            c2_PRS2_l_vl = self.interpLin(0.002, 0.005, 0.6, 0.1, paf)
                            c2_PRS3_l_vl = self.interpLin(0.003, 0.01, 0.6, 0.1, paf)
                            c2_PRS4_l_vl = self.interpLin(0.004, 0.012, 0.6, 0.1, paf)                            
                            c2_PRT2_g_vg = self.interpLin(0.003, 0.005, 0.6, 0.1, paf)
                            c2_PRT3_g_vg = self.interpLin(0.007, 0.045, 0.6, 0.1, paf)
                            c2_PRT4_g_vg = self.interpLin(0.008, 0.06, 0.6, 0.1, paf)
                            c2_PRT2_g_vl = self.interpLin(0.003, 0.005, 0.6, 0.1, paf)
                            c2_PRT3_g_vl = self.interpLin(0.009, 0.045, 0.6, 0.1, paf)
                            c2_PRT4_g_vl = self.interpLin(0.01, 0.06, 0.6, 0.1, paf)
                            c2_PRT2_l_vg = self.interpLin(0.000001, 0.004, 0.6, 0.1, paf)
                            c2_PRT3_l_vg = self.interpLin(0.000001, 0.005, 0.6, 0.1, paf)
                            c2_PRT4_l_vg = self.interpLin(0.000001, 0.006, 0.6, 0.1, paf)
                            c2_PRT2_l_vl = self.interpLin(0.002, 0.005, 0.6, 0.1, paf)
                            c2_PRT3_l_vl = self.interpLin(0.003, 0.01, 0.6, 0.1, paf)
                            c2_PRT4_l_vl = self.interpLin(0.004, 0.012, 0.6, 0.1, paf)
                            if rho_s >= 0.006:
                                c2_PRS2_vg = c2_PRS2_g_vg
                                c2_PRS3_vg = c2_PRS3_g_vg
                                c2_PRS4_vg = c2_PRS4_g_vg
                                c2_PRS2_vl = c2_PRS2_g_vl
                                c2_PRS3_vl = c2_PRS3_g_vl
                                c2_PRS4_vl = c2_PRS4_g_vl
                            elif rho_s <= 0.0005:
                                c2_PRS2_vg = c2_PRS2_l_vg
                                c2_PRS3_vg = c2_PRS3_l_vg
                                c2_PRS4_vg = c2_PRS4_l_vg
                                c2_PRS2_vl = c2_PRS2_l_vl
                                c2_PRS3_vl = c2_PRS3_l_vl
                                c2_PRS4_vl = c2_PRS4_l_vl
                            else:
                                c2_PRS2_vg = self.interpLin(c2_PRS2_l_vg, c2_PRS2_g_vg, 0.0005, 0.006, rho_s)
                                c2_PRS3_vg = self.interpLin(c2_PRS3_l_vg, c2_PRS3_g_vg, 0.0005, 0.006, rho_s)
                                c2_PRS4_vg = self.interpLin(c2_PRS4_l_vg, c2_PRS4_g_vg, 0.0005, 0.006, rho_s)   
                                c2_PRS2_vl = self.interpLin(c2_PRS2_l_vl, c2_PRS2_g_vl, 0.0005, 0.006, rho_s)
                                c2_PRS3_vl = self.interpLin(c2_PRS3_l_vl, c2_PRS3_g_vl, 0.0005, 0.006, rho_s)
                                c2_PRS4_vl = self.interpLin(c2_PRS4_l_vl, c2_PRS4_g_vl, 0.0005, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c2_PRT2_vg = c2_PRT2_g_vg
                                c2_PRT3_vg = c2_PRT3_g_vg
                                c2_PRT4_vg = c2_PRT4_g_vg
                                c2_PRT2_vl = c2_PRT2_g_vl
                                c2_PRT3_vl = c2_PRT3_g_vl
                                c2_PRT4_vl = c2_PRT4_g_vl
                            elif rho_t <= 0.0005:
                                c2_PRT2_vg = c2_PRT2_l_vg
                                c2_PRT3_vg = c2_PRT3_l_vg
                                c2_PRT4_vg = c2_PRT4_l_vg
                                c2_PRT2_vl = c2_PRT2_l_vl
                                c2_PRT3_vl = c2_PRT3_l_vl
                                c2_PRT4_vl = c2_PRT4_l_vl
                            else:
                                c2_PRT2_vg = self.interpLin(c2_PRT2_l_vg, c2_PRT2_g_vg, 0.0005, 0.006, rho_t)
                                c2_PRT3_vg = self.interpLin(c2_PRT3_l_vg, c2_PRT3_g_vg, 0.0005, 0.006, rho_t)
                                c2_PRT4_vg = self.interpLin(c2_PRT4_l_vg, c2_PRT4_g_vg, 0.0005, 0.006, rho_t)   
                                c2_PRT2_vl = self.interpLin(c2_PRT2_l_vl, c2_PRT2_g_vl, 0.0005, 0.006, rho_t)
                                c2_PRT3_vl = self.interpLin(c2_PRT3_l_vl, c2_PRT3_g_vl, 0.0005, 0.006, rho_t)
                                c2_PRT4_vl = self.interpLin(c2_PRT4_l_vl, c2_PRT4_g_vl, 0.0005, 0.006, rho_t)
                        if vC_s >= 0.5:
                            c2_PRS2 = c2_PRS2_vg
                            c2_PRS3 = c2_PRS3_vg
                            c2_PRS4 = c2_PRS4_vg
                        elif vC_s <= 0.25:
                            c2_PRS2 = c2_PRS2_vl
                            c2_PRS3 = c2_PRS3_vl
                            c2_PRS4 = c2_PRS4_vl
                        else:
                            c2_PRS2 = self.interpLin(c2_PRS2_vl, c2_PRS2_vg, 0.25, 0.5, vC_s)
                            c2_PRS3 = self.interpLin(c2_PRS3_vl, c2_PRS3_vg, 0.25, 0.5, vC_s)
                            c2_PRS4 = self.interpLin(c2_PRS4_vl, c2_PRS4_vg, 0.25, 0.5, vC_s)
                        if vC_t >= 0.5:
                            c2_PRT2 = c2_PRT2_vg
                            c2_PRT3 = c2_PRT3_vg
                            c2_PRT4 = c2_PRT4_vg
                        elif vC_t <= 0.25:
                            c2_PRT2 = c2_PRT2_vl
                            c2_PRT3 = c2_PRT3_vl
                            c2_PRT4 = c2_PRT4_vl
                        else:
                            c2_PRT2 = self.interpLin(c2_PRT2_vl, c2_PRT2_vg, 0.25, 0.5, vC_t)
                            c2_PRT3 = self.interpLin(c2_PRT3_vl, c2_PRT3_vg, 0.25, 0.5, vC_t)
                            c2_PRT4 = self.interpLin(c2_PRT4_vl, c2_PRT4_vg, 0.25, 0.5, vC_t)
                            
                    if any("3" in i for i in self.c_condition) or any("4" in i for i in self.c_condition):
                        c3_PRS1 = 0
                        c3_PRT1 = 0
                        if paf <= 0.1:
                            if rho_s >= 0.006:
                                c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 0, 0.000001, 0.045, 0.06
                            elif rho_s <= 0.0005:
                                c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 0, 0.000001, 0.005, 0.006
                            else:
                                c3_PRS2 = 0.000001
                                c3_PRS3 = self.interpLin(0.005, 0.045, 0.0005, 0.006, rho_s)
                                c3_PRS4 = self.interpLin(0.006, 0.06, 0.0005, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 0, 0.000001, 0.045, 0.06                               
                            elif rho_t <= 0.0005:
                                c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 0, 0.000001, 0.005, 0.006
                            else:
                                c3_PRT2 = 0.000001
                                c3_PRT3 = self.interpLin(0.005, 0.045, 0.0005, 0.006, rho_t)
                                c3_PRT4 = self.interpLin(0.006, 0.06, 0.0005, 0.006, rho_t)
                        if paf >= 0.6:
                            if rho_s >= 0.006:
                                c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 0, 0.000001, 0.007, 0.008
                            elif rho_s <= 0.0005:
                                c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 0, 0.000001, 0.000001, 0.000001
                            else:
                                c3_PRS2 = 0.000001
                                c3_PRS3 = self.interpLin(0.000001, 0.007, 0.0005, 0.006, rho_s)
                                c3_PRS4 = self.interpLin(0.000001, 0.008, 0.0005, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 0, 0.000001, 0.007, 0.008
                            elif rho_t <= 0.0005:
                                c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 0, 0.000001, 0.000001, 0.000001
                            else:
                                c3_PRT2 = 0.000001
                                c3_PRT3 = self.interpLin(0.000001, 0.007, 0.0005, 0.006, rho_t)
                                c3_PRT4 = self.interpLin(0.000001, 0.008, 0.0005, 0.006, rho_t)
                        else:
                            c3_PRS2 = 0.000001
                            c3_PRS3_g = self.interpLin(0.007, 0.045, 0.6, 0.1, paf)
                            c3_PRS4_g = self.interpLin(0.008, 0.06, 0.6, 0.1, paf)
                            c3_PRS3_l = self.interpLin(0.000001, 0.005, 0.6, 0.1, paf)
                            c3_PRS4_l = self.interpLin(0.000001, 0.006, 0.6, 0.1, paf)
                            c3_PRT2 = 0.000001
                            c3_PRT3_g = self.interpLin(0.007, 0.045, 0.6, 0.1, paf)
                            c3_PRT4_g = self.interpLin(0.008, 0.06, 0.6, 0.1, paf)
                            c3_PRT3_l = self.interpLin(0.000001, 0.005, 0.6, 0.1, paf)
                            c3_PRT4_l = self.interpLin(0.000001, 0.006, 0.6, 0.1, paf)
                            if rho_s >= 0.006:
                                c3_PRS3 = c3_PRS3_g
                                c3_PRS4 = c3_PRS4_g
                            elif rho_s <= 0.0005:
                                c3_PRS3 = c3_PRS3_l
                                c3_PRS4 = c3_PRS4_l
                            else:
                                c3_PRS3 = self.interpLin(c3_PRS3_l, c3_PRS3_g, 0.0005, 0.006, rho_s)
                                c3_PRS4 = self.interpLin(c3_PRS4_l, c3_PRS4_g, 0.0005, 0.006, rho_s)
                            if rho_t >= 0.006:
                                c3_PRT3 = c3_PRT3_g
                                c3_PRT4 = c3_PRT4_g
                            elif rho_t <= 0.0005:
                                c3_PRT3 = c3_PRT3_l
                                c3_PRT4 = c3_PRT4_l
                            else:
                                c3_PRT3 = self.interpLin(c3_PRT3_l, c3_PRT3_g, 0.0005, 0.006, rho_t)
                                c3_PRT4 = self.interpLin(c3_PRT4_l, c3_PRT4_g, 0.0005, 0.006, rho_t)
                        c4_PRS1 = c3_PRS1
                        c4_PRS2 = c3_PRS2
                        c4_PRS3 = c3_PRS3
                        c4_PRS4 = c3_PRS4
                        c4_PRT1 = c3_PRT1
                        c4_PRT2 = c3_PRT2
                        c4_PRT3 = c3_PRT3
                        c4_PRT4 = c3_PRT4
                        
                    # Specify dummy values for conditions that do not control
                    if not any("1" in i for i in self.c_condition):
                        c1_PRS1, c1_PRS2, c1_PRS3, c1_PRS4 = 1, 1, 1, 1
                        c1_PRT1, c1_PRT2, c1_PRT3, c1_PRT4 = 1, 1, 1, 1
                    if not any("2" in i for i in self.c_condition):
                        c2_PRS1, c2_PRS2, c2_PRS3, c2_PRS4 = 1, 1, 1, 1
                        c2_PRT1, c2_PRT2, c2_PRT3, c2_PRT4 = 1, 1, 1, 1
                    if not any("3" in i for i in self.c_condition):
                        c3_PRS1, c3_PRS2, c3_PRS3, c3_PRS4 = 1, 1, 1, 1
                        c3_PRT1, c3_PRT2, c3_PRT3, c3_PRT4 = 1, 1, 1, 1
                    if not any("4" in i for i in self.c_condition):
                        c4_PRS1, c4_PRS2, c4_PRS3, c4_PRS4 = 1, 1, 1, 1
                        c4_PRT1, c4_PRT2, c4_PRT3, c4_PRT4 = 1, 1, 1, 1

                    # Take the minimum outcomes among all specified conditions
                    PRS1, PRS2, PRS3, PRS4 = [
                        min(
                            c1_PRS1, c2_PRS1, c3_PRS1, c4_PRS1), min(
                            c1_PRS2, c2_PRS2, c3_PRS2, c4_PRS2), min(
                            c1_PRS3, c2_PRS3, c3_PRS3, c4_PRS3), min(
                            c1_PRS4, c2_PRS4, c3_PRS4, c4_PRS4)]
                    PRT1, PRT2, PRT3, PRT4 = [
                        min(
                            c1_PRT1, c2_PRT1, c3_PRT1, c4_PRT1), min(
                            c1_PRT2, c2_PRT2, c3_PRT2, c4_PRT2), min(
                            c1_PRT3, c2_PRT3, c3_PRT3, c4_PRT3), min(
                            c1_PRT4, c2_PRT4, c3_PRT4, c4_PRT4)]
                        
                else:
                    # Determine values required to calculate table constants 
                    # for ASCE 41-17 Tables 10-8 and 10-9 
                    # Note: Table 10-9 is the same for ASCE 41-17 and 41-13
                    
                    # Note e and general requirement
                    rho_t_s = min(max(self.Aswx / self.AgShearcx, 0.0005), 0.0175)
                    # Note e and general requirement
                    rho_t_t = min(max(self.Aswy / self.AgShearcy, 0.0005), 0.0175)
                    ss_ratio_s = max(0.2, self.ss_ratio_s)  # General requirement
                    ss_ratio_t = max(0.2, self.ss_ratio_t)  # General requirement
                    Naf_s = max((self.Ned * 1000) / (self.Ancx * self.fck),
                                0.1)  # Note b requirement
                    Naf_t = max((self.Ned * 1000) / (self.Ancy * self.fck),
                                0.1)  # Note b requirement
                    rho_l_s = self.Ascx / self.Agcx
                    rho_l_t = self.Ascy / self.Agcy
    
                    if any("3" in i for i in self.c_condition):
                        # Columns controlled by inadequate development or
                        # splicing along the clear height
    
                        # Compute constant "a"
                        if (self.shape == 'Rectangular Column') |\
                                (self.shape == 'Cruciform'):
                            a_s = max(
                                0.042 -
                                0.043 *
                                Naf_s +
                                0.63 *
                                rho_t_s -
                                0.023 *
                                ss_ratio_s,
                                0)  # Local s-direction
                            a_t = max(
                                0.042 -
                                0.043 *
                                Naf_t +
                                0.63 *
                                rho_t_t -
                                0.023 *
                                ss_ratio_t,
                                0)  # Local t-direction
                            b_s = max(a_s, 0.5 / (5 + Naf_s * self.fck /
                                                  (0.8 * rho_t_s * self.fywd)) - 0.01)
                            b_t = max(a_t, 0.5 / (5 + Naf_t * self.fck /
                                                  (0.8 * rho_t_t * self.fywd)) - 0.01)
                        elif self.shape == 'Circular':
                            a_s = max(
                                0.06 -
                                0.06 *
                                Naf_s +
                                1.3 *
                                rho_t_s -
                                0.037 *
                                ss_ratio_s,
                                0)  # Local s-direction
                            a_t = max(
                                0.06 -
                                0.06 *
                                Naf_t +
                                1.3 *
                                rho_t_t -
                                0.037 *
                                ss_ratio_t,
                                0)  # Local t-direction
                            b_s = max(a_s, 0.65 / (5 + Naf_s * self.fck /
                                                   (0.8 * rho_t_s * self.fywd)) - 0.01)
                            b_t = max(a_t, 0.65 / (5 + Naf_t * self.fck /
                                                   (0.8 * rho_t_t * self.fywd)) - 0.01)
    
                        if Naf_s > 0.5:  # Note a requirement
                            b_s5 = max(a_s, 0.5 / (5 + 0.5 * self.fck /
                                                   (0.8 * rho_t_s * self.fywd)) - 0.01)
                            b_s7 = max(a_s, 0.5 / (5 + 0.7 * self.fck /
                                                   (0.8 * rho_t_s * self.fywd)) - 0.01)
                            b_s = self.interpLin(b_s5, b_s7, 0.5, 0.7, Naf_s)
                        if Naf_t > 0.5:  # Note a requirement
                            b_t5 = max(a_t, 0.5 / (5 + 0.5 * self.fck /
                                                   (0.8 * rho_t_t * self.fywd)) - 0.01)
                            b_t7 = max(a_t, 0.5 / (5 + 0.7 * self.fck /
                                                   (0.8 * rho_t_t * self.fywd)) - 0.01)
                            b_t = self.interpLin(b_t5, b_t7, 0.5, 0.7, Naf_t)
    
                        # Compute constant "c"
                        c_s = max(0.24 - 0.4 * Naf_s, 0)
                        c_t = max(0.24 - 0.4 * Naf_t, 0)
    
                        # Determine rotation limits
                        PRS1, PRS2, PRS3, PRS4 = 0, min(
                            0.005, 0.15 * a_s), 0.5 * b_s, 0.7 * b_s
                        PRT1, PRT2, PRT3, PRT4 = 0, min(
                            0.005, 0.15 * a_t), 0.5 * b_t, 0.7 * b_t
    
                    else:
                        # Columns not controlled by inadequate development or
                        # splicing along the clear height
                        
                        # Compute constant "a"
                        a_s = min(max(0, (rho_t_s * self.fywd) /
                                      (8 * rho_l_s * self.fyk)), 0.025)
                        a_t = min(max(0, (rho_t_t * self.fywd) /
                                      (8 * rho_l_t * self.fyk)), 0.025)
                        # Compute constant "b"
                        # Note e requirement
                        b_s = min(max(0, 0.012 - 0.085 * Naf_s + 12 *
                                      min(rho_t_s, 0.0075), a_s), 0.06)
                        # Note e requirement
                        b_t = min(max(0, 0.012 - 0.085 * Naf_t + 12 *
                                      min(rho_t_t, 0.0075), a_t), 0.06)
    
                        # Compute constant "c"
                        c_s = min(0.15 + 36 * rho_t_s, 0.4)
                        c_t = min(0.15 + 36 * rho_t_t, 0.4)
    
                        # Determine rotation limits (must be in ascending order)
                        PRS1, PRS2, PRS3, PRS4 = 0, 0.1 * b_s, 0.5 * b_s, 0.7 * b_s
                        PRT1, PRT2, PRT3, PRT4 = 0, 0.1 * b_t, 0.5 * b_t, 0.7 * b_t

        else:
            PRS1, PRS2, PRS3, PRS4 = 0, 0, 0, 0
            PRT1, PRT2, PRT3, PRT4 = 0, 0, 0, 0

        PRS = [PRS1, PRS2, PRS3, PRS4]
        PRT = [PRT1, PRT2, PRT3, PRT4]
        PRS = [round(i, 6) for i in PRS]
        PRT = [round(i, 6) for i in PRT]

        return PRS
