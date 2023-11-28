# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 14:03:20 2017

@author: xinka
"""

import pandas as pd
import numpy as np
from gurobipy import *


def blend_optimize(price, quantity, product_spec, CC, confidence, No2 = False):
    # Raw materials
    raw_spec = pd.read_excel("Data\\Raw_spec_201901.xlsx", index_col=0)
    raw_spec['cost'] = price
    raw, High_Cu, Low_Cu, High_Zn, Low_Zn, High_Pb, Low_Pb, High_Sn, Low_Sn,High_Ni, Low_Ni,\
    High_Al, Low_Al, High_Mn, Low_Mn, High_Fe, Low_Fe, cost= multidict(raw_spec.T.to_dict('list'))
    
    # Product specs
    element, max_spec, min_spec = multidict({
        "Cu":[product_spec.High_Cu,product_spec.Low_Cu],
        "Zn":[product_spec.High_Zn,product_spec.Low_Zn],
        "Pb":[product_spec.High_Pb,product_spec.Low_Pb],
        "Sn":[product_spec.High_Sn,product_spec.Low_Sn],
        "Ni":[product_spec.High_Ni,product_spec.Low_Ni],
        "Al":[product_spec.High_Al,product_spec.Low_Al],
        "Mn":[product_spec.High_Mn,product_spec.Low_Mn],
        "Fe":[product_spec.High_Fe,product_spec.Low_Fe]       
    })
    
    # Production quanity
    prod = quantity
    
    # Wrapup confidence
    s = confidence*2 - 1
    
    # Model
    blend_model = Model('Blending Starter')
    
    # Decision variables
    raw_demand = blend_model.addVars(raw, name="raw_demand")
    
    # Objective function
    blend_model.setObjective(raw_demand.prod(cost), GRB.MINIMIZE)
    
    # Specs constraints: if CC, use chance-constrained; else use deterministic
    if CC:
        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Cu[i] + Low_Cu[i])/2 + (High_Cu[i] - Low_Cu[i])/2*s) for i in raw) <= max_spec['Cu'] * prod, "spec_Cu_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Cu[i] + Low_Cu[i])/2 - (High_Cu[i] - Low_Cu[i])/2*s) for i in raw) >= min_spec['Cu'] * prod, "spec_Cu_lo")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Zn[i] + Low_Zn[i])/2 + (High_Zn[i] - Low_Zn[i])/2*s) for i in raw) <= max_spec['Zn'] * prod, "spec_Zn_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Zn[i] + Low_Zn[i])/2 - (High_Zn[i] - Low_Zn[i])/2*s) for i in raw) >= min_spec['Zn'] * prod, "spec_Zn_lo")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Pb[i] + Low_Pb[i])/2 + (High_Pb[i] - Low_Pb[i])/2*s) for i in raw) <= max_spec['Pb'] * prod, "spec_Pb_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Pb[i] + Low_Pb[i])/2 - (High_Pb[i] - Low_Pb[i])/2*s) for i in raw) >= min_spec['Pb'] * prod, "spec_Pb_lo")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Sn[i] + Low_Sn[i])/2 + (High_Sn[i] - Low_Sn[i])/2*s) for i in raw) <= max_spec['Sn'] * prod, "spec_Sn_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Sn[i] + Low_Sn[i])/2 - (High_Sn[i] - Low_Sn[i])/2*s) for i in raw) >= min_spec['Sn'] * prod, "spec_Sn_lo")        

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Ni[i] + Low_Ni[i])/2 + (High_Ni[i] - Low_Ni[i])/2*s) for i in raw) <= max_spec['Ni'] * prod, "spec_Ni_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Ni[i] + Low_Ni[i])/2 - (High_Ni[i] - Low_Ni[i])/2*s) for i in raw) >= min_spec['Ni'] * prod, "spec_Ni_lo")        

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Al[i] + Low_Al[i])/2 + (High_Al[i] - Low_Al[i])/2*s) for i in raw) <= max_spec['Al'] * prod, "spec_Al_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Al[i] + Low_Al[i])/2 - (High_Al[i] - Low_Al[i])/2*s) for i in raw) >= min_spec['Al'] * prod, "spec_Al_lo")        

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Mn[i] + Low_Mn[i])/2 + (High_Mn[i] - Low_Mn[i])/2*s) for i in raw) <= max_spec['Mn'] * prod, "spec_Mn_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Mn[i] + Low_Mn[i])/2 - (High_Mn[i] - Low_Mn[i])/2*s) for i in raw) >= min_spec['Mn'] * prod, "spec_Mn_lo")        

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Fe[i] + Low_Fe[i])/2 + (High_Fe[i] - Low_Fe[i])/2*s) for i in raw) <= max_spec['Fe'] * prod, "spec_Fe_up")

        blend_model.addConstr(
                quicksum(raw_demand[i] * ((High_Fe[i] + Low_Fe[i])/2 - (High_Fe[i] - Low_Fe[i])/2*s) for i in raw) >= min_spec['Fe'] * prod, "spec_Fe_lo")        
        
        
    else:
        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Cu[i] + Low_Cu[i])/2 for i in raw), min_spec['Cu'] * prod, max_spec['Cu'] * prod, "spec_Cu")
    
        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Zn[i] + Low_Zn[i])/2 for i in raw), min_spec['Zn'] * prod, max_spec['Zn'] * prod, "spec_Zn")
    
        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Pb[i] + Low_Pb[i])/2 for i in raw), min_spec['Pb'] * prod, max_spec['Pb'] * prod, "spec_Pb")
    
        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Sn[i] + Low_Sn[i])/2 for i in raw), min_spec['Sn'] * prod, max_spec['Sn'] * prod, "spec_Sn")

        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Ni[i] + Low_Ni[i])/2 for i in raw), min_spec['Ni'] * prod, max_spec['Ni'] * prod, "spec_Ni")

        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Al[i] + Low_Al[i])/2 for i in raw), min_spec['Al'] * prod, max_spec['Al'] * prod, "spec_Al")

        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Mn[i] + Low_Mn[i])/2 for i in raw), min_spec['Mn'] * prod, max_spec['Mn'] * prod, "spec_Mn")

        blend_model.addRange(
                quicksum(raw_demand[i] * (High_Fe[i] + Low_Fe[i])/2 for i in raw), min_spec['Fe'] * prod, max_spec['Fe'] * prod, "spec_Fe")
    
    
    # Mass constraint
    blend_model.addConstr(
        quicksum(raw_demand[i] for i in raw) >= prod, "mass")
    
    
    if No2 == False:
        blend_model.addConstr(
                raw_demand['No.2'] == 0, "no #2")
    
    
    # Optimize!
    blend_model.update()
    blend_model.setParam( 'OutputFlag', False )
    blend_model.optimize()
    
    # Return results
    demand = blend_model.getAttr('x', raw_demand)
    return(demand)