# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 21:10:25 2021

@author: STA94720
"""

import pandas as pd
import numpy as np

def leapfrog_to_dyna(input_csv, output='model.key'):
    
    '''
    This function reads in data from LeapFrog and produces an LS-DYNA keyword
    file that contains PART, solid_data, and NODE data. Intended for
    porting 3D soil mesh data. 
    
    Parameters
    ----------     
    input_csv : str
        Full file path for the input LeapFrog .csv data file
        Only filename required if located in cwd
    output : str
        Full file path desired for the output file (.key)
        Only filename required if located in cwd
        Default = 'model.key'
    
    Returns
    -------
    {file}.key : file
        LS-DYNA keyword file 
    '''

    # read in data
    df = pd.read_csv(input_csv) 
    
    # shift model coordinates so one corner aligns with the origin
    minX = df['X'].min()
    minY = df['Y'].min()
    minZ = df['Z'].min()
    df['X'] = df['X']-minX
    df['Y'] = df['Y']-minY
    df['Z'] = df['Z']-minZ
    
    ## NODE ##
    node_data = df[['Id','X','Y','Z']]
    node_data['TC'] = int(0)
    node_data['RC'] = int(0)
    node_data_str = node_data.to_string(header=False, index=False, index_names=False).split('\n')
    NODE = [','.join(ele.split()) for ele in node_data_str]
    NODE = '\n'.join(NODE)
    
    ## PART ##
    geo_title = df.columns[-1].title()
    parts = pd.DataFrame(df[df.columns[-1]].unique())
    parts = parts.rename(columns = {parts.columns[-1]:geo_title})
    parts['PID'] = parts.index + 1
    df = df.merge(parts[[df.columns[-1].title(),'PID']],how='left')
    part_data = parts[['PID']]
    part_data['SID'] = 1
    part_data['MID'] = part_data['PID']
    part_data[['EOSID','HGID','GRAV','ADPOPT','TMID']] = 0
    part_data_str = [None] * len(part_data)
    for i in range(0,len(part_data)):
        part_data_str[i] = '*PART\n' +\
        str(parts[geo_title].loc[i]) + '\n' +\
        ','.join(pd.DataFrame(part_data.loc[i]).to_string(header=False, index=False, index_names=False).split('\n')) + '\n$'
    PART = '\n'.join(part_data_str)
    
    ## ELEMENT_SOLID ##
    # instantiate columns for other solid element nodes
    df['N2'] = df['Id']
    df['N3'] = df['Id']
    df['N4'] = df['Id']
    df['N5'] = df['Id']
    df['N6'] = df['Id']
    df['N7'] = df['Id']
    df['N8'] = df['Id']
    
    # populate solid element nodes with correct values assuming current Id is N1
    df_search = df.set_index(['X','Y','Z'])
    for i in range(0,len(df)):
        N2x = df['X'].loc[i] + df['dX'].loc[i]
        N2y = df['Y'].loc[i]
        N2z = df['Z'].loc[i]
        try:
            N2 = int(df_search.at[(N2x,N2y,N2z),'Id'])
        except KeyError:
            N2 = np.nan
        df['N2'].loc[i] = N2
        
        N3x = df['X'].loc[i] + df['dX'].loc[i]
        N3y = df['Y'].loc[i] + df['dY'].loc[i]
        N3z = df['Z'].loc[i]
        try:
            N3 = int(df_search.at[(N3x,N3y,N3z),'Id'])
        except KeyError:
            N3 = np.nan
        df['N3'].loc[i] = N3
        
        N4x = df['X'].loc[i]
        N4y = df['Y'].loc[i] + df['dY'].loc[i]
        N4z = df['Z'].loc[i]
        try:
            N4 = int(df_search.at[(N4x,N4y,N4z),'Id'])
        except KeyError:
            N4 = np.nan
        df['N4'].loc[i] = N4
        
        N5x = df['X'].loc[i]
        N5y = df['Y'].loc[i]
        N5z = df['Z'].loc[i] + df['dZ'].loc[i]
        try:
            N5 = int(df_search.at[(N5x,N5y,N5z),'Id'])
        except KeyError:
            N5 = np.nan
        df['N5'].loc[i] = N5
        
        N6x = df['X'].loc[i] + df['dX'].loc[i]
        N6y = df['Y'].loc[i]
        N6z = df['Z'].loc[i] + df['dZ'].loc[i]
        try:
            N6 = int(df_search.at[(N6x,N6y,N6z),'Id'])
        except KeyError:
            N6 = np.nan
        df['N6'].loc[i] = N6
        
        N7x = df['X'].loc[i] + df['dX'].loc[i]
        N7y = df['Y'].loc[i] + df['dY'].loc[i]
        N7z = df['Z'].loc[i] + df['dZ'].loc[i]
        try:
            N7 = int(df_search.at[(N7x,N7y,N7z),'Id'])
        except KeyError:
            N7 = np.nan
        df['N7'].loc[i] = N7
        
        N8x = df['X'].loc[i]
        N8y = df['Y'].loc[i] + df['dY'].loc[i]
        N8z = df['Z'].loc[i] + df['dZ'].loc[i]
        try:
            N8 = int(df_search.at[(N8x,N8y,N8z),'Id'])
        except KeyError:
            N8 = np.nan
        df['N8'].loc[i] = N8
        
    # create dataframe for solid_data data
    solid_data = pd.DataFrame(df.index + 1)
    solid_data = solid_data.rename(columns = {solid_data.columns[-1]:'EID'})
    solid_data['PID'] = df['PID']
    solid_data['N1'] = df['Id']
    solid_data['N2'] = df['N2']
    solid_data['N3'] = df['N3']
    solid_data['N4'] = df['N4']
    solid_data['N5'] = df['N5']
    solid_data['N6'] = df['N6']
    solid_data['N7'] = df['N7']
    solid_data['N8'] = df['N8']
    solid_data = solid_data.dropna()
    solid_data['EID'] = solid_data.index + 1
    solid_data = solid_data.astype(int)
    
    solid_data_str = solid_data.to_string(header=False, index=False, index_names=False).split('\n')
    ELEMENT_SOLID = [','.join(ele.split()) for ele in solid_data_str]
    ELEMENT_SOLID = '\n'.join(ELEMENT_SOLID)
    
    # write keyword file
    key_data = '*KEYWORD\n$\n'+\
        '*NODE\n' + NODE + '\n$\n' +\
        '*ELEMENT_SOLID\n' + ELEMENT_SOLID + '\n$\n' +\
        PART + '\n$\n' +\
        '*END'
    with open('model.key','w') as file:
        file.write(key_data)