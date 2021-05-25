# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 13:01:57 2018

@author: kevin.stanton

This class writes MAT_HYSTERETIC_BEAM card given output from class Shear_Calcs
"""

import math
from decimal import Decimal
from .util import util
import os


class Elastic_Beam():

    def __init__(self, aData):

        # Initialize section area
        self.area = aData.section_area

    def __process__(self, aData):

        # Shape of section
        shape = aData.section_type

        # Define moment of inertia
        if shape == 'Circular':
            # Diameter of circlular section (m)
            D = float(aData.diameter_circ) / 1000.0
            # Uncracked moment of inertia of the section (m4)
            I_xx = math.pi / 4 * (D / 2)**4
            I_yy = math.pi / 4 * (D / 2)**4
            # Cracked moment of inertia of the section (m4)
            I_xx_crack = I_xx / 3
            I_yy_crack = I_yy / 3

        elif (shape == 'Rectangular Beam') |\
             (shape == 'Rectangular Column'):
            # Get the width of the section (m)
            w = float(aData.rect_width) / 1000
            # Get the height of the section (m)
            h = float(aData.rect_height) / 1000
            # Uncracked moment of inertia of the section (m4)
            I_xx = w * h**3 / 12
            I_yy = h * w**3 / 12
            # Cracked moment of inertia of the section (m4)
            I_xx_crack = I_xx / 3
            I_yy_crack = I_yy / 3

        elif shape == 'Cruciform':
            # Get the widths of the section (m)
            bX = float(aData.bX) / 1000
            bY = float(aData.bX) / 1000
            # Get the heights of the section (m)
            hY = float(aData.tY) / 1000
            hX = float(aData.tX) / 1000
            # Uncracked moment of inertia of the section (m4)
            I_xx = bX * hY**3 / 12
            I_yy = bY * hX**3 / 12
            # Cracked moment of inertia of the section (m4)
            I_xx_crack = I_xx / 3
            I_yy_crack = I_yy / 3

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
            # Uncracked moment of inertia of the section (m4)
            I_xx = b * H * (y_cog - H / 2)**2 + b * H**3 / 12 + \
                h * B * (H + h / 2 - y_cog)**2 + h**3 * B / 12
            I_yy = b**3 * H / 12 + B**3 * h / 12
            # Cracked moment of inertia of the section (m4)
            I_xx_crack = I_xx / 3
            I_yy_crack = I_yy / 3

        # Define the modulus of elasticity of concrete (Pa)
        E = float(aData.Ec) * 1e6

        # Return output arguments
        return E, I_xx, I_yy, I_xx_crack, I_yy_crack

    def write_mat(self, aData):

        # Initialize curves
        E, I_xx, I_yy, I_xx_crack, I_yy_crack = self.__process__(aData)
        print('\n   Final card data:')
        print('         Young\'s modulus = ' + str(round(E, 0)))
        print('         I_xx = ' + str(round(I_xx, 6)))
        print('         I_xx (cracked) = ' + str(round(I_xx_crack, 6)))
        print('         I_yy = ' + str(round(I_yy, 6)))
        print('         I_yy (cracked) = ' + str(round(I_yy_crack, 6)))
        print('         Section area = ' + str(round(self.area, 5)) + '\n')

        # Input parameters
        zero_10 = '0'.rjust(10)

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
        iss = str(float(I_xx_crack))
        itt = str(float(I_yy_crack))
        area = ("{:.3E}".format(Decimal(area))).rjust(10)
        itt = ("{:.3E}".format(Decimal(itt))).rjust(10)
        iss = ("{:.3E}".format(Decimal(iss))).rjust(10)
        j = '0.0'.rjust(10)
        sa = '0.0'.rjust(10)
        ist = '0.0'.rjust(10)

        # Define material card parameters
        if not aData.pSave:
            keyData = []
            keyData.append('*KEYWORD\n')
            keyData.append('$\n')
            keyData.append('*MAT_ELASTIC_TITLE\n')
            keyData.append(util.getName(aData) + '\n')
            keyData.append(
                '%s%s%s%s%s%s%s%s\n' %
                (mid, ro, e, pr, zero_10, zero_10, zero_10, zero_10))
            keyData.append('$\n')

            # Write section card
            keyData.append('*SECTION_BEAM_TITLE\n')
            keyData.append(util.getName(aData) + '\n')
            keyData.append('%s%s%s%s%s%s%s\n' %
                           (sid, elform, shrf, qririd, cst, scoor, nsm))
            keyData.append('%s%s%s%s%s%s\n' % (area, iss, itt, j, sa, ist))
            keyData.append('$\n')
            keyData.append('*END\n')
        else:
            fwrite = open(
                aData.pSave +
                '\\' +
                util.getName(
                    aData,
                    '.key'),
                'w')
            fwrite.write('*KEYWORD\n')
            fwrite.write('$\n')
            fwrite.write('*MAT_ELASTIC_TITLE\n')
            fwrite.write(util.getName(aData) + '\n')
            fwrite.write('%s%s%s%s%s%s%s%s\n' %
                         (mid, ro, e, pr, zero_10, zero_10, zero_10, zero_10))
            fwrite.write('$\n')

            # Write section card
            fwrite.write('*SECTION_BEAM_TITLE\n')
            fwrite.write(util.getName(aData) + '\n')
            fwrite.write('%s%s%s%s%s%s%s\n' %
                         (sid, elform, shrf, qririd, cst, scoor, nsm))
            fwrite.write('%s%s%s%s%s%s\n' % (area, iss, itt, j, sa, ist))
            fwrite.write('$\n')
            fwrite.write('*END\n')

            # Close file handle
            fwrite.close()

        # Write output in JSON format
        if not aData.pSave == False:
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
        jDat['value'] = ''.join(keyData)

        return jDat
