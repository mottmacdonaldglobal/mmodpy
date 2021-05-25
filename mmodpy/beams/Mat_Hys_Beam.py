# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:01:57 2018

@author: kevin.stanton

This class writes MAT_HYSTERETIC_BEAM card given output from class Shear_Calcs
"""

import pandas as pd
import math
import numpy as np
from .Shear_Calcs import Shear_Calcs
from decimal import Decimal
from .biLinearize import biLin
from .util import util
import matplotlib.pyplot as plt
import os


class Mat_Hys_Beam(Shear_Calcs):

    # Envelope is a Boolean input
    def __init__(
            self,
            aData,
            section_area,
            envelope=False,
            n_strainlim=0,
            bilinearflag=False):
        print('         Envelope = ' + str(envelope))

        # Read data in a dataframe format
        # self.P_m must be replaced in order to do the enveloping method
        if not(envelope):
            self.P_m = pd.read_csv(aData.pSave + '\\' + str(aData.file_name) + '.txt',
                                   delim_whitespace=True,
                                   error_bad_lines=False,
                                   skiprows=15).iloc[:, [12, 14, 31, 34]]
            self.M_phi = pd.read_csv(aData.pSave + '\\' + str(aData.file_name) + '.txt',
                                     delim_whitespace=True,
                                     error_bad_lines=False,
                                     skiprows=15).iloc[:, [50, 53, 67, 69]].dropna()
            self.mPhi_long = pd.read_csv(aData.pSave + '\\' + str(aData.file_name) + '.txt',
                                         delim_whitespace=True,
                                         error_bad_lines=False,
                                         skiprows=15).iloc[:, [50, 53, 41, 43, 42, 67, 69, 57, 59, 58]].dropna()
        else:
            # Read from separate CSV
            self.P_m = pd.read_csv(aData.pSave +
                                   '\\' +
                                   str('Moment_Interaction_Envelope' +
                                       aData.file_name +
                                       '.csv'))
            colind = [12, 15, 29, 31]
            colstep = 38
            ind_mphi = [x + colstep * n_strainlim for x in colind]
            self.M_phi = pd.read_csv(aData.pSave + '\\' + str(aData.file_name) + '.txt',
                                     delim_whitespace=True,
                                     error_bad_lines=False,
                                     skiprows=15).iloc[:, ind_mphi].dropna()
            colind = [12, 15, 3, 5, 4, 29, 31, 19, 21, 20]
            colstep = 38
            ind_mphi = [x + colstep * n_strainlim for x in colind]
            self.mPhi_long = pd.read_csv(aData.pSave + '\\' + str(aData.file_name) + '.txt',
                                         delim_whitespace=True,
                                         error_bad_lines=False,
                                         skiprows=15).iloc[:, ind_mphi].dropna()

        # Add column names
        if envelope:
            self.P_m = self.P_m.reindex(columns=['Px', 'Mx', 'Py', 'My'])

        self.P_m.columns = ['N (N)', 'Mxx (N-m)', 'N (N)', 'Myy (N-m)']
        self.M_phi.columns = ['Mxx (N-m)', 'phi_xx (1/m)', 'Myy (N-m)',
                              'phi_yy (1/m)']
        self.mPhi_long.columns = [
            'Mxx (N-m)',
            'phi_xx (1/m)',
            'ConfinedMinStrain_xx',
            'SteelMinStrain_xx',
            'SteelMaxStrain_xx',
            'Myy (N-m)',
            'phi_yy (1/m)',
            'ConfinedMinStrain_yy',
            'SteelMinStrain_yy',
            'SteelMaxStrain_yy']

        # Convert floating points to integers
        self.P_m = self.P_m.astype(int)

        # Initialize section area
        self.area = section_area

        # Read in plot settings
        if aData.devMode == 'on':
            plot_fig = True
        else:
            plot_fig = False

        # Get bi-linearized and cleaned M_phi curves
        if aData.cLimits == 'explicit':
            self.xP, self.yP, self.M_ult, self.k_ult, ind = util.expYield(
                self.mPhi_long, aData)
            self.ind_xx = ind[0]  # x-curvature index at Yield
            self.ind_yy = ind[1]  # y-curvature index at Yield
        else:
            if aData.cLimits == 'strict':
                constrained = True
            else:
                constrained = False
            self.xP, self.yP, M_phi_bilin, self.M_ult, self.k_ult = biLin(
                self.M_phi, aData.file_name, plot_fig=plot_fig, constrained=constrained, cleanCurves=aData.cleanCurves)
            self.ind_xx = np.searchsorted(np.squeeze(
                self.M_phi[['phi_xx (1/m)']]), self.xP[0])[0]  # x-curvature index at Yield
            self.ind_yy = np.searchsorted(np.squeeze(
                self.M_phi[['phi_yy (1/m)']]), self.xP[1])[0]  # y-curvature index at Yield

        # Write moment curvature to csv file
        self.M_phi.to_csv(aData.pSave + '\\' + 'Moment_Curvature_RAW_' +
                          aData.file_name + '.csv', index=False)

        # Clean moment curvature plots according to settings
        xM, xPhi = self.M_phi[self.M_phi.columns[0]
                              ], self.M_phi[self.M_phi.columns[1]]
        yM, yPhi = self.M_phi[self.M_phi.columns[2]
                              ], self.M_phi[self.M_phi.columns[3]]
        xM_raw, xPhi_raw = xM, xPhi
        yM_raw, yPhi_raw = yM, yPhi
        if aData.cleanCurves:
            xPhi, xM, numCutX = util.cutOutlier(xPhi, xM)
            yPhi, yM, numCutY = util.cutOutlier(yPhi, yM)
            self.ind_xx = self.ind_xx + numCutX
            self.ind_yy = self.ind_yy + numCutY

        # Plot interpretation of My if necessary
        if plot_fig:
            if aData.cLimits == 'explicit':
                plt.figure()
                plt.title('M_phi_xx-' + str(aData.file_name))
                plt.plot(xPhi, xM, 'o')
                plt.plot(xPhi_raw, xM_raw, marker='*')
                plt.plot(self.xP[0], self.yP[0], color='red', marker='^')
                plt.savefig(aData.pSave + '\\' + 'M_phi_xx-' +
                            str(aData.file_name) + '.png')
                plt.figure()
                plt.title('M_phi_yy-' + str(aData.file_name))
                plt.plot(yPhi, yM, 'o')
                plt.plot(yPhi_raw, yM_raw, marker='*')
                plt.plot(self.xP[1], self.yP[1], color='red', marker='^')
                plt.savefig(aData.pSave + '\\' + 'M_phi_yy-' +
                            str(aData.file_name) + '.png')

        # Implement an option for the user to specify whether they want to use
        # the fully bilinearized M-phi curve or just the point closest to the yield point from the
        # bilinearzed M-phi curve and XTRACT M-phi curve for all points
        # thereafter ("direct" method)
        if bilinearflag:
            print('         M-phi interpretation = bilinear')
            self.M_phi = M_phi_bilin
        else:
            print('         M-phi interpretation = direct (' + str(aData.cLimits) + ')')
            Mphi_xx = pd.DataFrame({'Mxx (N-m)': xM, 'phi_xx (1/m)': xPhi})
            Mphi_yy = pd.DataFrame({'Mxx (N-m)': yM, 'phi_xx (1/m)': yPhi})
            if len(Mphi_xx) < len(Mphi_yy):
                self.M_phi = pd.concat(
                    [Mphi_xx, Mphi_yy], axis=1, join_axes=[Mphi_xx.index])
            else:
                self.M_phi = pd.concat(
                    [Mphi_xx, Mphi_yy], axis=1, join_axes=[Mphi_yy.index])

        # Write curves in ".csv" files
        self.P_m.to_csv(aData.pSave + '\\' + 'Moment_Interaction_' +
                        aData.file_name + '.csv', index=False)
        self.M_phi.to_csv(aData.pSave + '\\' + 'Moment_Curvature_' +
                          aData.file_name + '.csv', index=False)

    def __process_curves__(self, aData):

        # Axial load (N)
        print('         Axial load = ' + str(aData.N_load) + 'kN')
        N = float(aData.N_load) * 1000

        # Shape of section
        shape = aData.section_type

        # Plastic hinge length (m)
        pl_hi_xx, pl_hi_yy = util.getPlasticHingeLength(aData)

        # Define moment of inertia
        if shape == 'Circular':
            # Diameter of circlular section (m)
            D = float(aData.diameter_circ) / 1000.0
            # Get the uncracked moment of inertia of the section (m4)
            I_xx = math.pi / 4 * (D / 2)**4
            I_yy = math.pi / 4 * (D / 2)**4

        elif (shape == 'Rectangular Beam') |\
             (shape == 'Rectangular Column'):
            # Get the width of the section (m)
            w = float(aData.rect_width) / 1000
            # Get the height of the section (m)
            h = float(aData.rect_height) / 1000
            # Get the uncracked moment of inertia of the section (m4)
            I_xx = w * h**3 / 12
            I_yy = h * w**3 / 12

        elif shape == 'T Beam':
            # Get the width of the flange (m)
            B = float(aData.f_width) / 1000
            # Get the height of the flange (m)
            h = float(aData.f_thick) / 1000
            # Get the width of the web (m)
            b = float(aData.t_width) / 1000
            # Get the height of the web (m)
            H = float(aData.t_height) / 1000 - float(aData.f_thick) / 1000
            # Get the center of gravity (m)
            A = B * h + H * b
            y_cog = ((H + h / 2) * h * B + H**2 * b / 2) / A
            # Get the uncracked moment of inertia of the section (m4)
            I_xx = b * H * (y_cog - H / 2)**2 + b * H**3 / 12 + \
                h * B * (H + h / 2 - y_cog)**2 + h**3 * B / 12
            I_yy = b**3 * H / 12 + B**3 * h / 12

        elif shape == 'Cruciform':
            # Total section width (m)
            bX = float(aData.bX) / 1000
            # Total section height (m)
            bY = float(aData.bY) / 1000
            # Horizontal flange width (m)
            tX = float(aData.tX) / 1000
            # Vertical flange width (m)
            tY = float(aData.tX) / 1000
            # Eccentricities (m)
            cX = float(aData.cX) / 1000
            cY = float(aData.cY) / 1000
            # Get the uncracked moment of inertia of the section (m4)
            I_xx = bX**3 * tX / 12 + (bY - tX)**3 * tY / 12
            I_yy = bY**3 * tY / 12 + (bX - tY)**3 * tX / 12

        # Define uncracked section properties
        self.I_noCrack_xx = I_xx
        self.I_noCrack_yy = I_yy

        # Define curves for each direction
        P_m_xx = self.P_m.iloc[:, :2]
        P_m_yy = self.P_m.iloc[:, 2:]

        P_m_xx = P_m_xx.sort_values(by=['Mxx (N-m)'], ascending=False)
        P_m_yy = P_m_yy.sort_values(by=['Myy (N-m)'], ascending=False)

        M_phi_xx = self.M_phi.iloc[:, :2]
        M_phi_yy = self.M_phi.iloc[:, 2:]

        # Calculate the moment plastic rotation curve
        def __M_theta_pl__(P_m, M_phi, E, pl_hi, N, sectionname, axis):

            # Slice the moment interaction diagram (full to half)
            P_m = P_m.iloc[:int(len(P_m) - (len(P_m) - 1) / 2), :]

            # Find compressive capacity based on give M-N diagram
            Ncom = max(P_m.iloc[:, 0])
            print('         Compressive capacity = ' + str(Ncom / 1000) + 'kN')

            # Find tensile capacity based on give M-N diagram
            Nten = min(P_m.iloc[:, 0])
            print('         Tensile capacity = ' + str(Nten / 1000) + 'kN')

            # Take appropriate outputs from bi-linearization
            if axis == 'xx':
                cl_ind = self.ind_xx
                phi_y = self.xP[0]
                My = self.yP[0]
                M_ult = self.M_ult[0]
                k_ult = self.k_ult[0]
                EI_noCrack = self.I_noCrack_xx * E
                I_noCrack = self.I_noCrack_xx
            else:
                cl_ind = self.ind_yy
                phi_y = self.xP[1]
                My = self.yP[1]
                M_ult = self.M_ult[1]
                k_ult = self.k_ult[1]
                EI_noCrack = self.I_noCrack_yy * E
                I_noCrack = self.I_noCrack_yy

            # Cracked section properties
            EI = My / phi_y
            I = EI / E
            ratio = I / I_noCrack

            # Get plastic curvature from the M-phi curve - φpl (1/m)
            print('         I_concrete (uncracked) = ' + str(round(I_noCrack, 6)))
            print('         I_concrete+steel (cracked) = ' + str(round(I, 6)))
            print('         EI_concrete (uncracked) = ' + str(round(EI_noCrack, 0)))
            print('         EI_concrete+steel (cracked) = ' + str(round(EI, 0)))
            print('         Yield curvature = ' + str(round(phi_y, 4)))

            # phi_plastic = phi - M/EI where EI = My/phi_y at yield point
            phi_pl = M_phi.values[:, 1] - M_phi.values[:, 0] / EI

            # Combine arrays with bending moment and plastic curvature
            M_phi_pl = np.vstack(
                (M_phi.values[:, 0], phi_pl)).transpose()[1:, :]

            # Drop plastic rotation rows with index < cl_ind
            M_phi_pl = M_phi_pl[cl_ind:, [0, 1]]

            # Drop plastic rotation rows with negative values (if they exist)
            M_phi_pl = M_phi_pl[M_phi_pl[:, 1] > 0]
            print('         My at specified N = ' + str(round(My, 0)))

            # Insert point where section yields
            M_phi_pl = np.vstack((np.array([My, 0]).reshape(1, 2), M_phi_pl))

            # Translate curve to start at phi_y = 0
            M_phi_pl[:, 1] = M_phi_pl[:, 1] - M_phi_pl[0, 1]

            # Remove repeated end values in M-Phi curve
            if M_phi_pl[len(M_phi_pl) - 1,
                        1] == M_phi_pl[len(M_phi_pl) - 2,
                                       1]:
                if aData.devMode == 'on':
                    print(
                        '         Removing repeated end values from moment curvature data...')
                M_phi_pl = M_phi_pl[:len(M_phi_pl) - 2, :]

            # Get plastic rotation - (radians)
            print('         Plastic hinge length = ' +
                  str(round(pl_hi, 5)) + 'm')
            theta_pl = M_phi_pl[:, 1] * pl_hi

            # Combine arrays with bending moment and plastic rotation
            M_theta_pl = np.vstack((M_phi_pl[:, 0], theta_pl)).transpose()

            # Convert arrays to dataframes
            M_phi_pl = pd.DataFrame(data=M_phi_pl,
                                    columns=['Moment (N-m)',
                                             'Plastic Curvature (1/m)'])
            M_theta_pl = pd.DataFrame(data=M_theta_pl,
                                      columns=['Moment (N-m)',
                                               'Plastic Rotation (radians)'])

            # Return moment - plastic rotation curve as well as other
            # properties
            return M_theta_pl, My, Ncom, Nten, I, phi_y, M_ult

        # Calculate moment - plastic rotation curve
        print('\n   Moment-plastic rotation (xx direction):')
        M_theta_pl_xx, My_xx, Ncom_xx, Nten_xx,\
            I_xx, phi_y_xx, M_ult_xx = __M_theta_pl__(P_m_xx, M_phi_xx, float(
                aData.Ec) * 1e6, pl_hi_xx, N, aData.file_name, 'xx')
        print('\n   Moment-plastic rotation (yy direction):')
        M_theta_pl_yy, My_yy, Ncom_yy, Nten_yy,\
            I_yy, phi_y_yy, M_ult_yy = __M_theta_pl__(P_m_yy, M_phi_yy, float(
                aData.Ec) * 1e6, pl_hi_yy, N, aData.file_name, 'yy')

        # Define the compressive and tensile axial capacity (N)
        Ncom = max(Ncom_xx, Ncom_yy)
        Nten = max(Nten_xx, Nten_yy)

        # Define the modulus of elasticity of concrete (Pa)
        E = float(aData.Ec) * 1e6

        # Concatenate curves
        M_theta_pl = pd.concat([M_theta_pl_xx, M_theta_pl_yy],
                               axis=1, join_axes=[M_theta_pl_xx.index])
        M_theta_pl.columns = ['Mxx (N-m)', 'ThetaP_xx (radians)',
                              'Myy (N-m)', 'ThetaP_yy (radians)']

        # Write curves in ".csv" files
        M_theta_pl.to_csv(aData.pSave + '\\' + 'Moment_Plastic_Rotation_' +
                          aData.file_name + '.csv', index=False)

        # Obtain shear capacity - axial load curve
        # Define a list with values of axial load
        incr = (Ncom - Nten) / 200
        axial_range = [incr * i for i in range(201)]
        axial_range = [Nten + i for i in axial_range]

        # Define the shear capacity in both directions for range
        VRdx, VRdy = [], []
        shear = Shear_Calcs(aData)
        for i in range(len(axial_range)):
            VRdx_temp, VRdy_temp = shear.Shear_cap(Ned=axial_range[i] / 1000)
            VRdx.append(VRdx_temp * 1000)
            VRdy.append(VRdy_temp * 1000)

        # Create arrays with shear - axial force values
        VRd_xx = np.array([axial_range, VRdx])
        VRd_yy = np.array([axial_range, VRdy])
        VRd_xx = pd.DataFrame(data=VRd_xx.transpose(),
                              columns=['N (N)', 'V_xx (N)'])
        VRd_yy = pd.DataFrame(data=VRd_yy.transpose(),
                              columns=['N (N)', 'V_yy (N)'])

        # Normalize arrays by axial force
        axial_range_n = [i / Nten if i <
                         0 else i / (-Ncom) for i in axial_range]
        VRd_xx_n, VRd_yy_n = VRd_xx, VRd_yy
        VRd_xx_n['N (N)'], VRd_yy_n['N (N)'] = axial_range_n, axial_range_n

        # Obtain shear force under expected axial load
        VRdx_ex, VRdy_ex = shear.Shear_cap(Ned=float(aData.N_load))
        VRdx_ex, VRdy_ex = VRdx_ex * 1000, VRdy_ex * 1000

        # Obtain shear force under zero axial load
        VRdx_noN, VRdy_noN = shear.Shear_cap(Ned=0)
        VRdx_noN, VRdy_noN = VRdx_noN * 1000, VRdy_noN * 1000

        # Return output arguments
        return M_theta_pl_xx, M_theta_pl_yy, My_xx, My_yy,\
            Ncom, Nten, E, VRd_xx, VRd_yy, VRd_xx_n, VRd_yy_n,\
            VRdx_ex, VRdy_ex, I_xx, I_yy, M_ult_xx, M_ult_yy,\
            VRdx_noN, VRdy_noN

    def write_mat(self, aData, PRS=[0, 0, 0, 0], PRT=[0, 0, 0, 0]):

        # Initialize curves
        M_theta_pl_xx, M_theta_pl_yy, My_xx, My_yy,\
            Ncom, Nten, E, VRd_xx, VRd_yy, VRd_xx_n,\
            VRd_yy_n, VRdx_ex, VRdy_ex, I_xx, I_yy, M_ult_xx,\
            M_ult_yy, VRdx_noN, VRdy_noN = self.__process_curves__(aData)
        print('\n   Final card data:')
        print('         Young\'s modulus = ' + str(round(E, 0)))
        print('         I_xx = ' + str(round(I_xx, 6)))
        print('         I_yy = ' + str(round(I_yy, 6)))
        print('         My_xx = ' + str(round(My_xx, 6)))
        print('         M_ult_xx = ' + str(round(M_ult_xx, 6)))
        print('         My_yy = ' + str(round(My_yy, 6)))
        print('         M_ult_yy = ' + str(round(M_ult_yy, 6)))
        print('         Section area = ' + str(round(self.area, 5)))
        print('         Shear capacity in local s (N=0) = ' + str(round(VRdx_noN, 0)))
        print('         Shear capacity in local t (N=0) = ' + str(round(VRdy_noN, 0)))
        print('         PRS = ' + str(PRS))
        print('         PRT = ' + str(PRT))
        print('         SS = [0, 0.0002, 0.0005, 0.001]')
        print('         ST = [0, 0.0002, 0.0005, 0.001]' + '\n')

        # Input parameters
        zero_10 = '0'.rjust(10)
        one_10 = '1'.rjust(10)

        # Determine area ratio (<1 if in-slab T-Beam)
        if aData.section_type == 'T-Beam' and aData.tbOption == 'in-slab':
            aEff = aData.section_area - aData.f_width * \
                aData.f_thick  # effective shear area
            f = aEff / aData.section_area  # factor to reduce density written to LS-DYNA
        else:
            f = 1

        # Parameters for first row of material card
        # Define the material id
        mid = str(aData.beam_num).rjust(10)
        # Define mass density of concrete (RO)
        ro = str(f * aData.den).rjust(10)
        # Define Young's modulus of concrete (E)
        e = "{:.1e}".format(E).replace('+', '').rjust(10)
        # Define Poisson's ratio of concrete (PR)
        pr = str(aData.v).rjust(10)
        # Abcissa definition for LCAT and LCAC
        iax = '2'.rjust(10)
        # Yield surface type of interaction
        isurf = '0'.rjust(10)
        # Isotropic hardening type
        ihard = '4'.rjust(10)
        # Flag for input of FEMA thresholds
        ifema = '3'.rjust(10)

        # Parameters for second row of material card
        # Yield moment vs rotation about local s-axis curve is (LCPMS)
        lcpms = str(10000 + (aData.beam_num) * 10).rjust(10)
        # Scale factor for yield moment force in LCPMS (SFS)
        sfs = str(int(M_ult_xx)).rjust(10)
        # Yield moment vs rotation about local t-axis curve is (LCPMT)
        lcpmt = str(10001 + (aData.beam_num) * 10).rjust(10)
        # Scale factor for yield moment force in LCPMT (SFT)
        sft = str(int(M_ult_yy)).rjust(10)
        # Axial tensile force vs strain curve id (LCAT)
        lcat = str(10002 + (aData.beam_num) * 10).rjust(10)
        # Scale factor for tensile force in LCAT (SFAT)
        sfat = str(abs(int(Nten))).rjust(10)
        # Axial compressive force vs strain curve id (LCAC)
        lcac = str(10003 + (aData.beam_num) * 10).rjust(10)
        # Scale factor for tensile force in LCAC (SFAC)
        sfac = str(abs(int(Ncom))).rjust(10)

        # Parameters for third row of material card
        # Moment interaction diagram for the local s-axis (ALPHA)
        pms = str(10004 + (aData.beam_num) * 10).rjust(10)
        pmsi = str(-10004 - (aData.beam_num) * 10).rjust(10)
        # Moment interaction diagram for the local t-axis (BETA)
        pmt = str(10005 + (aData.beam_num) * 10).rjust(10)
        pmti = str(-10005 - (aData.beam_num) * 10).rjust(10)
        # Pinching factor for flexural hysteresis (PINM)
        pinm = '0'.rjust(10)
        # Pinching factor for shear hysteresis (PINS)
        pins = '0'.rjust(10)
        # Location of plastic hinge 1 from node 1 (HLOC1)
        hloc1 = '0'.rjust(10)
        # Location of plastic hinge 2 from node 2 (HLOC2)
        hloc2 = '0'.rjust(10)

        # Parameters for fourth row of material card
        # Yield shear force vs inelastic shear strain (LCSHS)
        lcshs = str(10006 + (aData.beam_num) * 10).rjust(10)
        # Scale factor of yield shear force in the s-direction (SFSHS)
        sfshs = str(10007 + (aData.beam_num) * 10).rjust(10)
        sfshsi = str(-10007 - (aData.beam_num) * 10).rjust(10)
        # Yield shear force vs inelastic shear strain (LCSHS)
        lcsht = str(10008 + (aData.beam_num) * 10).rjust(10)
        # Scale factor of yield shear force in the t-direction (SFSHT)
        sfsht = str(10009 + (aData.beam_num) * 10).rjust(10)
        sfshti = str(-10009 - (aData.beam_num) * 10).rjust(10)
        # Define shear strength for expected axial force (N)
        vrdx = str(round(VRdx_ex, 0)).rjust(10)
        vrdy = str(round(VRdy_ex, 0)).rjust(10)
        # Define shear strength for 0 axial force
        vrdx0 = str(round(VRdx_noN, 0)).rjust(10)
        vrdy0 = str(round(VRdy_noN, 0)).rjust(10)

        # Parameters for eighth row of material card
        PRS = [str(round(i, 5)).rjust(10) for i in PRS]
        PRT = [str(round(i, 5)).rjust(10) for i in PRT]
        # Rotation limits in the local s-direction
        prs1 = PRS[0]
        prs2 = PRS[1]
        prs3 = PRS[2]
        prs4 = PRS[3]
        # Rotation limits in the local t-direction
        prt1 = PRT[0]
        prt2 = PRT[1]
        prt3 = PRT[2]
        prt4 = PRT[3]

        # Parameters for tenth row of material card
        # Shear strain thresholds in the local s-direction
        ss1 = zero_10
        ss2 = str(round(0.0002, 5)).rjust(10)
        ss3 = str(round(0.0005, 5)).rjust(10)
        ss4 = str(round(0.001, 5)).rjust(10)
        # Shear strain thresholds in the local t-direction
        st1 = zero_10
        st2 = str(round(0.0002, 5)).rjust(10)
        st3 = str(round(0.0005, 5)).rjust(10)
        st4 = str(round(0.001, 5)).rjust(10)

        # Define material card parameters
        fwrite = open(
            aData.pSave +
            '\\' +
            util.getName(
                aData,
                '.key'),
            'w')
        fwrite.write('*KEYWORD\n')
        fwrite.write('$\n')
        fwrite.write('*MAT_HYSTERETIC_BEAM_TITLE\n')
        fwrite.write(util.getName(aData) + '\n')
        fwrite.write('%s%s%s%s%s%s%s%s\n' %
                     (mid, ro, e, pr, iax, isurf, ihard, ifema))
        fwrite.write('%s%s%s%s%s%s%s%s\n' %
                     (lcpms, sfs, lcpmt, sft, lcat, sfat, lcac, sfac))
        fwrite.write('%s%s%s%s%s%s%s%s\n' %
                     (pmsi, pmti, zero_10, zero_10, pinm, pins, hloc1, hloc2))
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (zero_10,
             zero_10,
             zero_10,
             zero_10,
             lcshs,
             sfshsi,
             lcsht,
             sfshti))
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10))
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10))
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10))
        fwrite.write('%s%s%s%s%s%s%s%s\n' %
                     (prs1, prs2, prs3, prs4, prt1, prt2, prt3, prt4))
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10))
        fwrite.write('%s%s%s%s%s%s%s%s\n' %
                     (ss1, ss2, ss3, ss4, st1, st2, st3, st4))
        fwrite.write('$\n')

        # Define section card parameters
        # Card 1
        sid = mid  # same as material id
        elform = '2'.rjust(10)
        shrf = '0.0'.rjust(10)
        qririd = '0.0'.rjust(10)
        cst = '0.0'.rjust(10)
        scoor = '0.0'.rjust(10)
        nsm = '0.0'.rjust(10)

        # Card 2
        area = str(float(self.area))
        iss = str(float(I_xx))
        itt = str(float(I_yy))
        area = ("{:.3E}".format(Decimal(area))).rjust(10)
        itt = ("{:.3E}".format(Decimal(itt))).rjust(10)
        iss = ("{:.3E}".format(Decimal(iss))).rjust(10)
        j = '0.0'.rjust(10)
        sa = '0.0'.rjust(10)
        ist = '0.0'.rjust(10)

        # Write section card
        fwrite.write('*SECTION_BEAM_TITLE\n')
        fwrite.write(util.getName(aData) + '\n')
        fwrite.write('%s%s%s%s%s%s%s\n' %
                     (sid, elform, shrf, qririd, cst, scoor, nsm))
        fwrite.write('%s%s%s%s%s%s\n' % (area, iss, itt, j, sa, ist))
        fwrite.write('$\n')

        # Define curves
        # Define plastic moment vs plastic rotation curve (s-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Yield Moment vs Rotation about local s-axis (LCPMS)\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (lcpms,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        for i in range(0, len(M_theta_pl_xx)):
            fwrite.write('%s%s\n' % (str('%.4E' % Decimal(str(M_theta_pl_xx.iloc[i, 1]))).rjust(
                20), str('%.4E' % Decimal(str(M_theta_pl_xx.iloc[i, 0] / M_ult_xx))).rjust(20)))
        fwrite.write('$\n')

        # Define plastic moment vs plastic rotation curve (t-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Yield Moment vs Rotation about local t-axis (LCPMT)\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (lcpmt,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))

        for i in range(0, len(M_theta_pl_yy)):
            fwrite.write('%s%s\n' % (str('%.4E' % Decimal(str(M_theta_pl_yy.iloc[i, 1]))).rjust(
                20), str('%.4E' % Decimal(str(M_theta_pl_yy.iloc[i, 0] / M_ult_yy))).rjust(20)))

        fwrite.write('$\n')

        # Define axial tensile force vs strain curve (s-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Axial Tensile Force vs Strain (LCAT)\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (lcat,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        fwrite.write('                 0.0                 1.0\n')
        fwrite.write('                 0.1                 1.0\n')
        fwrite.write('$\n')

        # Define compressive tensile force vs strain curve (t-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Axial Compressive Force vs Strain (LCAC)\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (lcac,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        fwrite.write('                 0.0                 1.0\n')
        fwrite.write('                 0.1                 1.0\n')
        fwrite.write('$\n')

        # Define moment interaction curve (s-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Moment Interaction curve about local s-axis (ALPHA)\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (pms,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        # in DYNA, (+) Tension, (-) Compression.
        # in XTRACT, (+) Compression, (-) Tension
        # Multiply y axis of generated PM curve by -1
        for i in range(len(self.P_m)):
            fwrite.write('%s%s\n' % (str('%.4E' % Decimal(str(self.P_m.iloc[i, 1]))).rjust(
                20), str('%.4E' % Decimal(str(self.P_m.iloc[i, 0] * -1.0))).rjust(20)))
        fwrite.write('$\n')

        # Define moment interaction curve (t-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Moment Interaction curve about local t-axis (BETA)\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (pmt,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        # in DYNA, (+) Tension, (-) Compression.
        # in XTRACT, (+) Compression, (-) Tension
        # Multiply y axis of generated PM curve by -1
        for i in range(len(self.P_m)):
            fwrite.write('%s%s\n' % (str('%.4E' % Decimal(str(self.P_m.iloc[i, 3]))).rjust(
                20), str('%.4E' % Decimal(str(self.P_m.iloc[i, 2] * -1.0))).rjust(20)))
        fwrite.write('$\n')

        # Define yield shear vs inelastic shear strain curve (s-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write(
            'Yield Shear Force vs Inelastic Strain about local s-axis\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (lcshs,
             zero_10,
             zero_10,
             vrdx0,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        fwrite.write('                 0.0                 1.0\n')
        fwrite.write('                 0.1                 1.0\n')
        fwrite.write('$\n')

        # Define yield shear vs inelastic shear strain curve (t-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write(
            'Yield Shear Force vs Inelastic Strain about local t-axis\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (lcsht,
             zero_10,
             zero_10,
             vrdy0,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        fwrite.write('                 0.0                 1.0\n')
        fwrite.write('                 0.1                 1.0\n')
        fwrite.write('$\n')

        # Define shear vs axial force curve (s-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Shear Force vs Axial Force about local s-axis\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (sfshs,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        # The shear force is divided by the shear capacity VRdx_noN (vrdx0 string form)
        # Meaning unity at axial force = 0
        for i in reversed(range(len(VRd_xx_n))):
            fwrite.write('%s%s\n' % (str('%.4E' % Decimal(str(VRd_xx_n.iloc[i, 0]))).rjust(
                20), str('%.4E' % Decimal(str(VRd_xx_n.iloc[i, 1] / VRdx_noN))).rjust(20)))
        fwrite.write('$\n')

        # Define shear vs axial force curve (t-axis)
        fwrite.write('*DEFINE_CURVE_TITLE\n')
        fwrite.write('Shear Force vs Axial Force about local t-axis\n')
        fwrite.write(
            '%s%s%s%s%s%s%s%s\n' %
            (sfsht,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             zero_10,
             one_10,
             zero_10))
        # The shear force is divided by the shear capacity VRdy_noN (vrdy0 string form)
        # Meaning unity at axial force = 0
        for i in reversed(range(len(VRd_yy_n))):
            fwrite.write('%s%s\n' % (str('%.4E' % Decimal(str(VRd_yy_n.iloc[i, 0]))).rjust(
                20), str('%.4E' % Decimal(str(VRd_yy_n.iloc[i, 1] / VRdy_noN))).rjust(20)))
        fwrite.write('$\n')
        fwrite.write('*END\n')

        # Close file handle
        fwrite.close()

        # Write output in JSON format
        keyFile = open(
            os.path.join(
                aData.pSave,
                util.getName(
                    aData,
                    '.key')),
            'r')
        keyData = keyFile.read()
        keyFile.close()
        jDat = {}
        jDat['file_name'] = util.getName(aData, '.key')
        jDat['value'] = keyData

        return jDat
