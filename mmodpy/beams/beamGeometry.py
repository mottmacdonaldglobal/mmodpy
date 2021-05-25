# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 12:54:15 2018

@author: kevin.stanton

This class initializes an XTRACT model given user defined data from xtract3.py
"""


class beamGeometry():

    def __sect__(self, aData):

        # Write section builder properties for circular section
        if aData.section_type == 'Circular':
            section_area = 3.14 * pow(float(aData.diameter_circ) / 2000, 2.0)

        # Write section builder properties for rectangular beam section
        elif aData.section_type == 'Rectangular Beam':
            section_area = float(aData.rect_width) / 1000 * \
                float(aData.rect_height) / 1000

        # Write section builder properties for rectangular column section
        elif aData.section_type == 'Rectangular Column':
            section_area = round(
                (float(
                    aData.rect_width) /
                    1000) *
                float(
                    aData.rect_height) /
                1000,
                5)

        # Write section builder properties for T Beam section
        elif aData.section_type == 'T Beam':
            section_area = (float(aData.f_width) / 1000) * float(aData.f_thick) / 1000 + float(
                aData.t_width) / 1000 * (float(aData.t_height) / 1000 - float(aData.f_thick) / 1000)

        # Write section builder properties for Cruciform section
        elif aData.section_type == 'Cruciform':
            section_area = round(((float(aData.tX) *
                                   float(aData.bY) +
                                   float(aData.tY) *
                                   float(aData.bX) -
                                   float(aData.tX) *
                                   float(aData.tY)) /
                                  (1000**2)), 5)
        return section_area

    def runElastic(self, aData):

        # Calculate required geometric data
        section_area = self.__sect__(aData)

        return section_area
