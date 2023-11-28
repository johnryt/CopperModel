# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 09:56:09 2019

@author: Guest User
"""

import pandas as pd
import numpy as np
idx = pd.IndexSlice


def raw_specs_middle():
    
    # Importing raw specs and calculate middle
    raw_specs = pd.read_excel('Data\\Raw_spec_201901.xlsx', index_col=0)
    ele_index = [f.split('_')[1] for f in list(raw_specs.columns)][0::2]
    raw_med_specs = pd.DataFrame(0, index=raw_specs.index, columns=ele_index)
    
    for raw in raw_med_specs.index:
        high = raw_specs.loc[raw,:][0::2]
        low = raw_specs.loc[raw,:][1::2]
        med = (high.values+low.values)/2
        raw_med_specs.loc[raw,:] = med
    
    for raw in raw_med_specs.index:
        tot_ele = raw_med_specs.sum(axis=1).loc[raw]
        if tot_ele > 100:
            raw_med_specs.loc[raw, :] = raw_med_specs.loc[raw, :].div(tot_ele).mul(100)

    return(raw_med_specs)


def raw_price_historical_sim(s):
    # Importing raw historical prices and deflate
    raw_price_hist = pd.read_excel('Data\\Historical prices.xlsx', index_col=0)
    deflator = pd.read_excel('Data\\Deflator.xls', index_col=0, header=0)
    raw_med_specs = raw_specs_middle()
    
    for date in raw_price_hist.index:
        y = date.year
        deflate = deflator.loc[pd.datetime(y,1,1)]
        raw_price_hist.loc[date] = raw_price_hist.loc[date]/deflate.values*100
    
    ele_price_hist = raw_price_hist.iloc[:, 0:8]

    # Calculate prices back based on price formation model
    other_metal_count = s  # How much other metal's values are counted
    ele_price_hist_nominal = ele_price_hist.copy()
    ele_price_hist_nominal.iloc[:,2:] = ele_price_hist.iloc[:,2:].mul(other_metal_count)
    
    scrap_comps = raw_med_specs.iloc[10:,]
    scrap_metal_values = pd.DataFrame(np.matmul(scrap_comps.div(100), ele_price_hist_nominal.transpose()), 
                                      index=scrap_comps.index, columns=ele_price_hist_nominal.index).transpose()
        
    sp_low = 911 # Nominal value - consumer buy, Based on other_metal_count = 0
    sp_high = 1256 # Based on other_metal_count = 0
    alloy_price_init = scrap_metal_values.iloc[-1, :] - (sp_low + other_metal_count*(sp_high-sp_low))
    
    alloy_sim_index = raw_price_hist.iloc[-1,:].isna()
    raw_price_hist.loc[pd.datetime(2018,1,1), alloy_sim_index] = alloy_price_init.loc[alloy_sim_index]
    
    alloy_spread_diff_hist = ele_price_hist.loc[:, 'Ref_Cu'].diff().mul(0.424).fillna(0)
    alloy_spread_diff_hist_cum = np.cumsum(alloy_spread_diff_hist[::-1])[::-1].shift(-1).fillna(0)

    for s in alloy_sim_index[alloy_sim_index].index:
        alloy_spread = pd.Series(0, index=raw_price_hist.index)
        alloy_spread.iloc[-1] = scrap_metal_values.loc[pd.datetime(2018,1,1),s] - raw_price_hist.loc[pd.datetime(2018,1,1), s]
        alloy_spread.iloc[:] = alloy_spread.iloc[-1] - alloy_spread_diff_hist_cum
        raw_price_hist.loc[:, s] = scrap_metal_values.loc[:, s] - alloy_spread
        
    return (raw_price_hist)
    

def raw_price_future_sim(s, seed, vol_mean, vol_std, end_time=pd.datetime(2048,1,1)):
    raw_med_specs = raw_specs_middle()
    raw_price_future = pd.DataFrame(0, index=pd.date_range('20180101', end_time, freq='AS'), columns=raw_med_specs.index)

    # Initilize based on historical function
    price_2018 = raw_price_historical_sim(s).iloc[-1,:]
    raw_price_future.iloc[0,:] = price_2018
    ele_price_future = raw_price_future.iloc[:, 0:8]
    
    # Generate random walks
    np.random.seed(seed)
    vols = pd.DataFrame(np.random.normal(loc=vol_mean, scale=vol_std, size=(raw_price_future.shape[0], 8)), index=raw_price_future.index, columns=ele_price_future.columns)
    
    for i in range(vols.shape[0]-1):
        this_year_prices = ele_price_future.iloc[i, :]
        next_year_prices = this_year_prices.mul(vols.iloc[i+1, :]+1)
        ele_price_future.iloc[i+1, :] = next_year_prices
        
    raw_price_future.iloc[:, 0:8] = ele_price_future
    
    # Unalloyed scrap
    no1_spread_diff = raw_price_future.loc[:, 'Ref_Cu'].diff().mul(0.0648).fillna(0)
    no2_spread_diff = raw_price_future.loc[:, 'Ref_Cu'].diff().mul(0.185).fillna(0)
    
    no1_spread = pd.Series(0, index=raw_price_future.index, name='No.1')
    no1_spread.iloc[0] = raw_price_future.loc[pd.datetime(2018,1,1),'Ref_Cu'] - raw_price_future.loc[pd.datetime(2018,1,1),'No.1']
    no1_spread.iloc[:] = no1_spread.iloc[0] + no1_spread_diff.cumsum()
    raw_price_future.loc[:, 'No.1'] = raw_price_future.loc[:, 'Ref_Cu'] - no1_spread
    
    no2_spread = pd.Series(0, index=raw_price_future.index, name='No.2')
    no2_spread.iloc[0] = raw_price_future.loc[pd.datetime(2018,1,1),'Ref_Cu'] - raw_price_future.loc[pd.datetime(2018,1,1),'No.2']
    no2_spread.iloc[:] = no2_spread.iloc[0] + no2_spread_diff.cumsum()
    raw_price_future.loc[:, 'No.2'] = raw_price_future.loc[:, 'Ref_Cu'] - no2_spread
    
    # Alloyed scrap
    ele_price_future_nominal = ele_price_future.copy()
    ele_price_future_nominal.iloc[:,2:] = ele_price_future.iloc[:,2:].mul(s)
    
    scrap_comps = raw_med_specs.iloc[10:,]
    scrap_metal_values = pd.DataFrame(np.matmul(scrap_comps.div(100), ele_price_future_nominal.transpose()), 
                                      index=scrap_comps.index, columns=ele_price_future_nominal.index).transpose()
    
    alloy_spread_diff = raw_price_future.loc[:, 'Ref_Cu'].diff().mul(0.424).fillna(0)
    
    for s in scrap_metal_values.columns:
        alloy_spread = pd.Series(0, index=raw_price_future.index)
        alloy_spread.iloc[0] = scrap_metal_values.loc[pd.datetime(2018,1,1),s] - raw_price_future.loc[pd.datetime(2018,1,1), s]
        alloy_spread.iloc[:] = alloy_spread.iloc[0] + alloy_spread_diff.cumsum()
        raw_price_future.loc[:, s] = scrap_metal_values.loc[:, s] - alloy_spread
    
    return (raw_price_future)