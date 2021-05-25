# -*- coding: utf-8 -*-
"""
Created on Tue May 15 12:23:40 2018

@author: kevin.stanton
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


def biLin(
        M_phi,
        sectionName,
        plot_fig=False,
        constrained=False,
        cleanCurves=False):
    '''
    This method assigns a 2-stage piecewise linear curve fit to 2D array data.

    Inputs:
        M_phi = DataFrame of moment curvature values for x and y directions
        a, b = Guesses for slope (a) and intercept (b) on segment P1 (after yield)
        c, d = Guesses for slope (c) and intercept (d) on segment P2 (before yield)
    '''

    def pSegs(x, a, b, c, d):
        '''
        Function to compute piecewise segments.
        '''

        P1 = a * np.array(x) + b
        P2 = c * np.array(x) + d

        return np.minimum(P1, P2)

    def cutOutlier(x, y):

        i = 0
        xMod = []
        yMod = []
        while i < len(y) - 1:
            if float(y[i]) < float(y[i + 1]):
                xMod.append(x[i])
                yMod.append(y[i])
            i = i + 1

        return xMod, yMod

    def preFit(x, y):

        i = 1
        while i < len(y):
            if float(y[i]) == max(y):
                c = i
            i = i + 1

        x = x[:c]
        y = y[:c]

        return x, y

    def refineFit(x, y, index):

        i = 1
        c = len(y)
        while i <= len(y):
            if float(y[i]) > index:
                c = min(i + 4, len(y))
                i = len(y)
            i = i + 1

        x = x[:c]
        y = y[:c]

        return x, y

    ######################################
    ########### Bi-linearize: x ##########
    ######################################

    # Read in data
    y, x = M_phi.loc[:, 'Mxx (N-m)'], M_phi.loc[:, 'phi_xx (1/m)']

    # Remove outliers and prepare for the fitting process
    x, y = cutOutlier(x, y)
    xMod, yMod = preFit(x, y)

    # Set constants for initaial bi-linear guess
    a = max((yMod[len(yMod) - 1] - yMod[len(yMod) - 2]) /
            max(xMod[len(xMod) - 1] - xMod[len(xMod) - 2], 0.0001), 10)
    b = max(yMod) * 0.8
    c = (yMod[2] - yMod[1]) / (xMod[2] - xMod[1])
    d = 0
    if c < 0 or math.isnan(c):
        c = b / (max(x) / 4)

    # Fit piecewise function to identify the preliminary yield point
    pw0 = (a, b, c, d)
    pw, cov = curve_fit(pSegs, xMod, yMod, pw0)

    # Define preliminary yield point
    xP_xx = (pw[3] - pw[1]) / (pw[0] - pw[2])
    yP_xx = min(max(yMod), abs(
        pw[0] * ((pw[3] - pw[1]) / (pw[0] - pw[2])) + pw[1]))

    # Refine interpreation of the yield point
    a = (max(yMod) - yP_xx) / (max(xMod) - xP_xx)
    b = pw[1]
    c = yP_xx / xP_xx
    xMod, yMod = refineFit(xMod, yMod, yP_xx)
    sigma = np.ones(len(xMod))
    pw0 = (a, b, c, d)
    if constrained:
        sigma[[0, -1]] = 0.01
    pw, cov = curve_fit(pSegs, xMod, yMod, pw0, sigma=sigma)
    xP_xx = (pw[3] - pw[1]) / (pw[0] - pw[2])
    yP_xx = pw[0] * ((pw[3] - pw[1]) / (pw[0] - pw[2])) + pw[1]

    # Determine ultimate moment capacity coordinates
    M_ult_xx = max(y)
    k_ult_xx = np.interp(
        M_ult_xx, M_phi.loc[:, 'Mxx (N-m)'], M_phi.loc[:, 'phi_xx (1/m)'])

    # Convert result back to DF
    M_phi_biLin_xx = pd.DataFrame(
        {'Mxx (N-m)': pSegs(x, *pw), 'phi_xx (1/m)': x})

    # Remove non-monotonically increasing values
    M_phi_biLin_mI_xx = np.maximum.accumulate(M_phi_biLin_xx.values)

    M_phi_biLin_xx = pd.DataFrame(
        {'Mxx (N-m)': M_phi_biLin_mI_xx[:, 0], 'phi_xx (1/m)': M_phi_biLin_mI_xx[:, 1]})

    M_phi_biLin_xx = M_phi_biLin_xx.drop_duplicates()

    # =========================================================================
    # Plot results for sanity check
    if not cleanCurves:
        y, x = M_phi.loc[:, 'Mxx (N-m)'], M_phi.loc[:, 'phi_xx (1/m)']

    if plot_fig:
        plt.figure(1)
        plt.title('M_phi_xx')
        plt.plot(x, y, 'o', x, pSegs(x, *pw), '-')
        plt.plot(xP_xx, yP_xx, color='red', marker='^')
        plt.savefig('M_phi_xx-' + str(sectionName) + '.png')
        plt.show()

    # =========================================================================

    ######################################
    ########### Bi-linearize: y ##########
    ######################################
    y, x = M_phi.loc[:, 'Myy (N-m)'], M_phi.loc[:, 'phi_yy (1/m)']
    x, y = cutOutlier(x, y)
    xMod, yMod = preFit(x, y)

    # Set constants for initaial bi-linear guess
    a = max((yMod[len(yMod) - 1] - yMod[len(yMod) - 2]) /
            max(xMod[len(xMod) - 1] - xMod[len(xMod) - 2], 0.0001), 10)
    b = max(yMod) * 0.8
    c = (yMod[2] - yMod[0]) / (xMod[2] - xMod[0])
    if c < 0 or math.isnan(c):
        c = b / (max(x) / 4)

    # Fit piecewise function to identify the preliminary yield point
    pw0 = (a, b, c, d)
    pw, cov = curve_fit(pSegs, xMod, yMod, pw0)

    # Define preliminary yield point
    xP_yy = (pw[3] - pw[1]) / (pw[0] - pw[2])
    yP_yy = min(max(yMod), abs(
        pw[0] * ((pw[3] - pw[1]) / (pw[0] - pw[2])) + pw[1]))

    # Refine interpreation of the yield point
    a = (max(yMod) - yP_yy) / (max(xMod) - xP_yy)
    b = pw[1]
    c = yP_yy / xP_yy
    xMod, yMod = refineFit(xMod, yMod, yP_yy)
    sigma = np.ones(len(xMod))
    pw0 = (a, b, c, d)
    if constrained:
        sigma[[0, -1]] = 0.01
    pw, cov = curve_fit(pSegs, xMod, yMod, pw0, sigma=sigma)
    xP_yy = (pw[3] - pw[1]) / (pw[0] - pw[2])
    yP_yy = pw[0] * ((pw[3] - pw[1]) / (pw[0] - pw[2])) + pw[1]

    # Determine ultimate moment capacity coordinates
    M_ult_yy = max(y)
    k_ult_yy = np.interp(
        M_ult_yy, M_phi.loc[:, 'Myy (N-m)'], M_phi.loc[:, 'phi_yy (1/m)'])

    # Convert result back to DF
    M_phi_biLin_yy = pd.DataFrame(
        {'Myy (N-m)': pSegs(x, *pw), 'phi_yy (1/m)': x})

    # Remove non-monotonically increasing values
    M_phi_biLin_mI_yy = np.maximum.accumulate(M_phi_biLin_yy.values)
    M_phi_biLin_yy = pd.DataFrame(
        {'Myy (N-m)': M_phi_biLin_mI_yy[:, 0], 'phi_yy (1/m)': M_phi_biLin_mI_yy[:, 1]})
    M_phi_biLin_yy = M_phi_biLin_yy.drop_duplicates()

    # =========================================================================
    # Plot results for sanity check
    if not cleanCurves:
        y, x = M_phi.loc[:, 'Myy (N-m)'], M_phi.loc[:, 'phi_yy (1/m)']

    if plot_fig:
        plt.figure(2)
        plt.title('M_phi_yy')
        plt.plot(x, y, 'o', x, pSegs(x, *pw), '-')
        plt.plot(xP_yy, yP_yy, color='red', marker='^')
        plt.savefig('M_phi_yy-' + str(sectionName) + '.png')
        plt.show()

    # =========================================================================

    #########################################
    ###### Combine results for x and y ######
    #########################################

    if len(M_phi_biLin_xx.dropna()) < len(M_phi_biLin_yy.dropna()):
        M_phi_biLin = pd.concat([M_phi_biLin_xx, M_phi_biLin_yy], axis=1, join_axes=[
                                M_phi_biLin_xx.index])
    else:
        M_phi_biLin = pd.concat([M_phi_biLin_xx, M_phi_biLin_yy], axis=1, join_axes=[
                                M_phi_biLin_yy.index])

    xP = np.array([abs(xP_xx), abs(xP_yy)])  # phi
    yP = np.array([abs(yP_xx), abs(yP_yy)])  # M
    M_ult = np.array([abs(M_ult_xx), abs(M_ult_yy)])
    k_ult = np.array([abs(k_ult_xx), abs(k_ult_yy)])

    return xP, yP, M_phi_biLin, M_ult, k_ult
