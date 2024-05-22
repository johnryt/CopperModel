# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 09:47:53 2019

@author: Guest User
"""

import pandas as pd
from blending import blend_optimize
from raw_price_gen import raw_price_historical_sim, raw_price_future_sim


def raw_consumption(
        alloy_prod,
        alloy_breakdown,
        start_time = datetime(1999,1,1),
        end_time = datetime(2018,1,1),
        other_metal_count=1, seed=1120, vol_mean=0, vol_std=0.2, verbose=0):
    
    # Read specs
    prod_spec_with_breakdown = alloy_breakdown.copy()
    prod_spec_with_breakdown = prod_spec_with_breakdown.set_index('Primary code')
    prod_spec = prod_spec_with_breakdown.iloc[:, 2:18] # Subject to change
    prod_breakdown = prod_spec_with_breakdown.iloc[:, [1,-1]] # Subject to change
    raw_spec = pd.read_excel("Data\\Raw_spec_201901.xlsx", index_col=0)
    
    
    
    # Simulate prices
    raw_hist = raw_price_historical_sim(s=other_metal_count)
    raw_future = raw_price_future_sim(seed=seed, s=other_metal_count, vol_mean=vol_mean, vol_std=vol_std, end_time=end_time)
    raw_price = raw_hist.append(raw_future.iloc[1:,:])
    
    # Category separation according to ICA data
    pss = prod_spec.loc[prod_breakdown['Category'] == 'PSS']
    tube = prod_spec.loc[prod_breakdown['Category'] == 'Tube']
    rbs = prod_spec.loc[prod_breakdown['Category'] == 'RBS']
    wire = prod_spec.loc[prod_breakdown['Category'] == 'Wire']
    cast = prod_spec.loc[prod_breakdown['Category'] == 'Cast']
    
    total_raw_demand_series = pd.DataFrame(0.0, index=pd.date_range(start_time, end_time, freq='AS'), columns=raw_spec.index)
    
    for y in pd.date_range(start_time, end_time, freq='AS'):
        price = raw_price.loc[y,:]
        total_raw_demand = pd.Series(0, index=raw_price.columns)
    
        for i in range(len(tube)):
            product_spec = tube.iloc[i,:]
            alloy_name = tube.iloc[i,:].name
            tube_quantity = alloy_prod.loc[y,'tube']*prod_breakdown.loc[prod_breakdown['Category']=='Tube'].loc[alloy_name, 'Weight']
            raw_demand = blend_optimize(price, tube_quantity, product_spec, CC = True, confidence = 0.95)
            total_raw_demand += raw_demand.values()
        
        for i in range(len(rbs)):
            product_spec = rbs.iloc[i,:]
            alloy_name = rbs.iloc[i,:].name
            rbs_quantity = alloy_prod.loc[y,'RBS']*prod_breakdown.loc[prod_breakdown['Category']=='RBS'].loc[alloy_name, 'Weight']
            raw_demand = blend_optimize(price, rbs_quantity, product_spec, CC = True, confidence = 0.95)
            total_raw_demand += raw_demand.values()
        
        for i in range(len(pss)):
            product_spec = pss.iloc[i,:]
            alloy_name = pss.iloc[i,:].name
            pss_quantity = alloy_prod.loc[y,'PSS']*prod_breakdown.loc[prod_breakdown['Category']=='PSS'].loc[alloy_name, 'Weight']
            raw_demand = blend_optimize(price, pss_quantity, product_spec, CC = True, confidence = 0.95)
            total_raw_demand += raw_demand.values()
            
        for i in range(len(wire)):
            alloy_name = wire.iloc[i,:].name
            wire_quantity = alloy_prod.loc[y,'wire']*prod_breakdown.loc[prod_breakdown['Category']=='Wire'].loc[alloy_name, 'Weight']
            product_spec = wire.iloc[i,:]
            raw_demand = blend_optimize(price, wire_quantity, product_spec, CC = True, confidence = 0.95)
            total_raw_demand += raw_demand.values()
        
        for i in range(len(cast)):
            alloy_name = cast.iloc[i,:].name
            cast_quantity = alloy_prod.loc[y,'cast']*prod_breakdown.loc[prod_breakdown['Category']=='Cast'].loc[alloy_name, 'Weight']
            product_spec = cast.iloc[i,:]
            raw_demand = blend_optimize(price, cast_quantity, product_spec, CC = True, confidence = 0.95)
            total_raw_demand += raw_demand.values()
        
        total_raw_demand_series.loc[y,] = total_raw_demand
        if verbose == 1:
            print('Year', y.year, 'Optimized')
            
    return(total_raw_demand_series)