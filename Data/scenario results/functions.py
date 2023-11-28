import brewer2mpl as b2m
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
bmap = b2m.get_map('Dark2','qualitative',5)
# bmap = b2m.get_map('Paired','qualitative',4)
colors = bmap.mpl_colors
# from pylab import *
import matplotlib as mpl
# b2m.print_maps()
color_08 = 'darkblue'
color_10 = 'tab:green'
from datetime import datetime
from dateutil.relativedelta import relativedelta
idx = pd.IndexSlice

c1 = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']
c2 = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ffffb3']
c3 = ['#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e','#003c30']
c4 = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2']
c5 = ['#d73027','#f46d43','#fdae61','#fee090', '#abd9e9','#74add1','#4575b4']
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=c5)

def init_plot(fontsize=16,linewidth=3,figwidth=8,ax=0,figheight=5.5):
    params = {
       'axes.labelsize': fontsize,
       'font.size': fontsize,
       'legend.fontsize': fontsize,
       'xtick.labelsize': fontsize,
       'ytick.labelsize': fontsize,
       'text.usetex': False,
       'figure.figsize': [figwidth, figheight],
       'lines.linewidth': linewidth,
       'legend.framealpha': 1,
       'legend.frameon': False,
       'mathtext.default': 'regular'
       }

    mpl.rcParams.update(params)
    mpl.rcParams['axes.spines.left'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.bottom'] = False
    mpl.rcParams['axes.axisbelow'] = True
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['axes.grid.axis'] = 'y'
    mpl.rcParams['grid.color'] = '0.9'
    mpl.rcParams['grid.linewidth'] = 1
    mpl.rcParams['grid.linestyle'] = '-'
#     %config InlineBackend.figure_format ='retina'

def pplot(y, label, color='tab:pink', x=0, fontsize=16, linewidth=3, firstcall=0):
    if firstcall==1:
        plt.figure()
        params = {
           'axes.labelsize': fontsize,
           'font.size': fontsize,
           'legend.fontsize': fontsize,
           'xtick.labelsize': fontsize,
           'ytick.labelsize': fontsize,
           'text.usetex': False,
           'figure.figsize': [8, 5.5],
           'lines.linewidth': linewidth
           }
        mpl.rcParams.update(params)
        plt.rc('axes', axisbelow=True)
        plt.axes(frameon=0)
        plt.grid(axis='y',color='0.9',linestyle='-',linewidth=1)

    plt.plot(y, color = color, label = label, linewidth = linewidth)
    legend = plt.legend()
    frame = legend.get_frame()
    frame.set_facecolor('1.0')
    frame.set_edgecolor('1.0')
#     %config InlineBackend.figure_format ='retina'
    
def psubplot(ax, y, label, color='tab:pink', x=0, fontsize=16, linewidth=3, firstcall=0, width = 16):
    '''Takes the axes, pandas series, and label and gives prettier subplot, width recommended in multiples of 8 for number of subplots'''
    if firstcall == 1:
        params = {
        'axes.labelsize': 12,
        'font.size': 12,
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'text.usetex': False,
        'figure.figsize': [width, 5.5],
        'lines.linewidth': 3
        }
        mpl.rcParams.update(params)
        mpl.rcParams['axes.spines.left'] = False
        mpl.rcParams['axes.spines.right'] = False
        mpl.rcParams['axes.spines.top'] = False
        mpl.rcParams['axes.spines.bottom'] = False
#     %config InlineBackend.figure_format ='retina'

    ax.plot(y, label=label,color=color,linewidth=linewidth)
    legend = ax.legend()
    frame = legend.get_frame()
    frame.set_facecolor('1.0')
    frame.set_edgecolor('1.0')
    ax.grid(axis='y',color='0.9',linestyle='-',linewidth=1)
    ax.set_axisbelow(True)
    

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
    from mpl_toolkits.axes_grid1 import AxesGrid

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def panel_aic(model):
    k=1+len(model.params)
    L=abs(model.loglik)
    aic = 2*k-2*np.log(L)
    return aic

def hausman(fe, re):
    """
    Compute hausman test for fixed effects/random effects models
    b = beta_fe
    B = beta_re
    From theory we have that b is always consistent, but B is consistent
    under the alternative hypothesis and efficient under the null.
    The test statistic is computed as
    z = (b - B)' [V_b - v_B^{-1}](b - B)
    The statistic is distributed z \sim \chi^2(k), where k is the number
    of regressors in the model.
    Parameters
    ==========
    fe : statsmodels.regression.linear_panel.PanelLMWithinResults
        The results obtained by using sm.PanelLM with the
        method='within' option.
    re : statsmodels.regression.linear_panel.PanelLMRandomResults
        The results obtained by using sm.PanelLM with the
        method='swar' option.
    Returns
    =======
    chi2 : float
        The test statistic
    df : int
        The number of degrees of freedom for the distribution of the
        test statistic
    pval : float
        The p-value associated with the null hypothesis
    Notes
    =====
    The null hypothesis supports the claim that the random effects
    estimator is "better". If we reject this hypothesis it is the same
    as saying we should be using fixed effects because there are
    systematic differences in the coefficients.
    """
    
    import numpy.linalg as la
    
    # Pull data out
    b = fe.params
    B = re.params
    v_b = fe.cov
    v_B = re.cov

    # NOTE: find df. fe should toss time-invariant variables, but it
    #       doesn't. It does return garbage so we use that to filter
    df = b[np.abs(b) < 1e8].size

    # compute test statistic and associated p-value
    chi2 = np.dot((b - B).T, la.inv(v_b - v_B).dot(b - B))
    pval = stats.chi2.sf(chi2, df)

    return chi2, df, pval

def displacement_estimate(sd_scenario, scrap_baseline, 
                          mining_baseline):
    scrap_supply_scenario=sd_scenario.loc['19600101':, 'Scrap production']
    mining_supply_scenario=sd_scenario.loc['19600101':, ['Concentrate production', 'SX-EW production']].sum(axis=1)
    shock_scenario=scrap_supply_scenario-scrap_baseline
    mining_response_scenario=mining_supply_scenario-mining_baseline
    displacement_scenario=mining_response_scenario.cumsum().div(shock_scenario.cumsum()).mul(-1)
    
    dis_results=pd.DataFrame({'Shock': shock_scenario, 'Mining response': mining_response_scenario, 
                              'Displacement': displacement_scenario, 'Scrap %Increase': shock_scenario.cumsum()/scrap_baseline.cumsum(),
                              'Mining %Increase': mining_response_scenario.cumsum()/mining_baseline.cumsum(),
                              'Scrap %Increase 2020': shock_scenario.cumsum()/scrap_baseline.loc['20200101'],
                              'Mining %Increase 2020': mining_response_scenario.cumsum()/mining_baseline.loc['20200101']})
    return dis_results

def calculate_impacts(sd_simulated_cn, sd_simulated_rw, sd_simulated, mine_life_stats_all_tp,
                     direct_melt_demand_cn, direct_melt_demand_rw):
    '''Calculates environmental impacts based on production, returns total impacts, unit impacts, and quantities'''
    idx = pd.IndexSlice
    mine_life_stats_all = mine_life_stats_all_tp.transpose()
    sxew_id_operating = [i for i in mine_life_stats_all.columns.levels[0] if 
                    round(mine_life_stats_all.loc[:, idx[i,'Recovered metal production (kt)']].mean()/
                     mine_life_stats_all.loc[:, idx[i,'Paid metal production (kt)']].mean(),6) <= 1
                     and 'Inc_' not in str(i)]
    conc_id_operating = [i for i in mine_life_stats_all.columns.levels[0] if 
                    i not in sxew_id_operating and 'Inc_' not in str(i)]
    sxew_id_new = [i for i in mine_life_stats_all.columns.levels[0] if 
                    round(mine_life_stats_all.loc[:, idx[i,'Recovered metal production (kt)']].mean()/
                     mine_life_stats_all.loc[:, idx[i,'Paid metal production (kt)']].mean(),6) <= 1
                     and 'Inc_' in str(i)]
    conc_id_new = [i for i in mine_life_stats_all.columns.levels[0] if 
                    i not in sxew_id_new and 'Inc_' in str(i)]
    id_new = [i for i in mine_life_stats_all.columns.levels[0] if 'Inc_' in str(i)]
    id_op = [i for i in mine_life_stats_all.columns.levels[0] if 'Inc_' not in str(i)]
    mine_life_stats_panel_operating = mine_life_stats_all.loc[:,idx[id_op,:]].copy()
    mine_life_stats_panel_new = mine_life_stats_all.loc[:,idx[id_new,:]].copy()
    
    # Loading data:
    import os
    directory = os.path.dirname(os.path.dirname(os.getcwd()))
    semis_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Semis Distribution', index_col=0).loc[:'Latin America',:'Fraction Semis Production']
    semis_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Semis', index_col=0)
    mod_semis_distribution = semis_distribution / (1 - semis_distribution.loc['China'])

    refinery_primary_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Primary Refinery Ecoinvent', index_col=0).iloc[:,1:]
    refinery_secondary_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Scrap Ecoinvent', index_col=0, skiprows = 12).loc[:,:'Water use (m3)']
    refinery_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Primary Refinery Ecoinvent', index_col=0).iloc[:,0]
    mod_refinery_distribution = refinery_distribution / (1 - refinery_distribution.loc['China'])

    conc_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    conc_mult = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Mult', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    conc_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).iloc[:,0]
    conc_scale = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).loc[:,'TRACI Scale':]
    sxew_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Ecoinvent', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    sxew_mult = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Mult', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    sxew_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Ecoinvent', index_col=0).iloc[:,0]/100
    sxew_scale = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Ecoinvent', index_col=0).loc[:,'TRACI Scale':]
    mining_evolution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).loc[:,['Fraction','Identified Resources','Undiscovered Resources']]
    mine_eqn_coefficients = pd.DataFrame([[1.5777, -0.626], [2.0621, -1.208], [15.697, -0.573], [36.529, -0.351], 
                                          [0.073378, -0.094], [0.041232, -0.34]], 
                                         index = ['Conc TRACI', 'SX-EW TRACI','Conc Energy','SX-EW Energy',
                                                  'Conc Water','SX-EW Water'], columns = ['A','B'])
    mine_eqn_coefficients.index.name = '(A*ore_grade^B) * scale'

    scrap_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Scrap Ecoinvent', index_col=0, skiprows = 1).loc[:'No1',:'Water use (m3)']
    # direct_melt_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Semis Distribution', index_col=0).loc[:'Latin America',:'Fraction Semis Production']
    refined_metal_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Ref_Metals', index_col=0).iloc[:,2:]
    refined_metal_impacts.index = ['Ref_Zn','Ref_Pb','Ref_Sn','Ref_Ni','Ref_Al','Ref_Mn','Ref_Fe']
    mine_regionality = pd.read_excel(directory+'/Data/primary supply/Operating mine pool - countries.xlsx', index_col=0)

    # Mining evolution
    new_conc_dist = pd.DataFrame(np.array((list(mining_evolution.loc[:,'Fraction'])*23)).reshape(23, 7), index = mine_life_stats_panel_operating.index, columns = mining_evolution.index)
    new_sxew_dist = pd.DataFrame(np.array((list(sxew_distribution)*23)).reshape(23, 7), index = mine_life_stats_panel_operating.index, columns = sxew_distribution.index)

    # Other initialization
    # Semis first (excluding refined metals, those will be last so we have that year's copper impacts)
    countries = refinery_primary_impacts.index
    impacts = refinery_primary_impacts.columns
    traci_impacts = impacts[:10]
    energy_impacts = impacts[10:-1]
    water_impacts = list([impacts[-1]])
    colss = pd.MultiIndex.from_product([countries, impacts])
    idxs = pd.MultiIndex.from_product([direct_melt_demand_cn.index, \
                ['Direct Melt Scrap', 'Direct Melt Other', 'Fabrication', 'New SX-EW Mines',
                 'New Traditional Mines', 'Operating SX-EW Mines', 'Operating Traditional Mines', 
                 'Primary Refinery', 'Secondary Refinery', 'Refined Metals']])
    country_impacts = pd.DataFrame(0, index = idxs, columns = colss)
    country_unit_impacts = pd.DataFrame(0, index = idxs, columns = colss)
    country_quantities = pd.DataFrame(0, index = idxs, columns = countries)

    # Mine and SX-EW preliminary calculations
    op_mines_traci = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    new_mines_traci = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    op_sxew_traci = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    new_sxew_traci = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    op_mines_energy = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    new_mines_energy = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    op_sxew_energy = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    new_sxew_energy = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    op_mines_water = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    new_mines_water = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    op_sxew_water = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)
    new_sxew_water = pd.DataFrame(0, index = direct_melt_demand_cn.loc['20180101':].index, columns = countries)

    for c in countries:
        ####### Mining #######
        conc_op_region_id = [i for i in conc_id_operating if c in mine_regionality.loc[i,'Region']]
        sxew_op_region_id = [i for i in sxew_id_operating if c in mine_regionality.loc[i,'Region']]
        conc_cutoff_traci = 10
        sxew_cutoff_traci = 28
        conc_cutoff_energy = 150
        sxew_cutoff_energy = 450
        conc_cutoff_water = 0.27
        sxew_cutoff_water = 0.1

        # Operating Concentrate Mines - TRACI
        intermediate1 = (mine_eqn_coefficients.loc['Conc TRACI','A'] * mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc TRACI','B']) * conc_scale.loc[c,'TRACI Scale']
        intermediate1[intermediate1>conc_cutoff_traci] = conc_cutoff_traci
        op_mines_traci.loc[:,c] = (intermediate1 * \
         mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New Concentrate Mines - TRACI
        intermediate2 = ((mine_eqn_coefficients.loc['Conc TRACI','A'] * mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc TRACI','B']) * conc_scale.loc[c,'TRACI Scale']).apply(\
         lambda x: x * new_conc_dist.loc[:,c])
        intermediate2[intermediate2>conc_cutoff_traci] = conc_cutoff_traci
        new_mines_traci.loc[:,c] = (intermediate2 * \
         mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating SX-EW - TRACI
        intermediate3 = (mine_eqn_coefficients.loc['SX-EW TRACI','A'] * mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW TRACI','B']) * sxew_scale.loc[c,'TRACI Scale']
        intermediate3[intermediate3>sxew_cutoff_traci] = sxew_cutoff_traci
        op_sxew_traci.loc[:,c] = (intermediate3 * \
         mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New SX-EW - TRACI
        intermediate4 = ((mine_eqn_coefficients.loc['SX-EW TRACI','A'] * mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW TRACI','B']) * sxew_scale.loc[c,'TRACI Scale']).apply(\
         lambda x: x * new_sxew_dist.loc[:,c])
        intermediate4[intermediate4>sxew_cutoff_traci] = sxew_cutoff_traci
        new_sxew_traci.loc[:,c] = (intermediate4 * \
         mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating Concentrate Mines - Energy
        intermediate5 = (mine_eqn_coefficients.loc['Conc Energy','A'] * mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Energy','B'] ) * conc_scale.loc[c,'Energy Scale']
        intermediate5[intermediate5>conc_cutoff_energy] = conc_cutoff_energy
        op_mines_energy.loc[:,c] = (intermediate5 * \
         mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New Concentrate Mines - Energy
        intermediate6 = ((mine_eqn_coefficients.loc['Conc Energy','A'] * mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Energy','B']) * conc_scale.loc[c,'Energy Scale']).apply(\
         lambda x: x * new_conc_dist.loc[:,c])
        intermediate6[intermediate6>conc_cutoff_energy] = conc_cutoff_energy
        new_mines_energy.loc[:,c] = (intermediate6 * \
         mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating SX-EW - Energy
        intermediate7 = (mine_eqn_coefficients.loc['SX-EW Energy','A'] * mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Energy','B']) * sxew_scale.loc[c,'Energy Scale']
        intermediate7[intermediate7>sxew_cutoff_energy] = sxew_cutoff_energy
        op_sxew_energy.loc[:,c] = (intermediate7 * \
         mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New SX-EW - Energy
        intermediate8 = ((mine_eqn_coefficients.loc['SX-EW Energy','A'] * mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Energy','B']) * sxew_scale.loc[c,'Energy Scale']).apply(\
         lambda x: x * new_sxew_dist.loc[:,c])
        intermediate8[intermediate8>sxew_cutoff_energy] = sxew_cutoff_energy
        new_sxew_energy.loc[:,c] = (intermediate8 * \
         mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating Concentrate Mines - Water
        intermediate9 = (mine_eqn_coefficients.loc['Conc Water','A'] * mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Water','B'] ) * conc_scale.loc[c,'Water Scale']
        intermediate9[intermediate9>conc_cutoff_water] = conc_cutoff_water
        op_mines_water.loc[:,c] = (intermediate9 * \
         mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New Concentrate Mines - Water
        intermediate10 = ((mine_eqn_coefficients.loc['Conc Water','A'] * mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Water','B']) * conc_scale.loc[c,'Water Scale']).apply(\
         lambda x: x * new_conc_dist.loc[:,c])
        intermediate10[intermediate10>conc_cutoff_water] = conc_cutoff_water
        new_mines_water.loc[:,c] = (intermediate10 * \
         mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating SX-EW - Water
        intermediate11 = (mine_eqn_coefficients.loc['SX-EW Water','A'] * mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Water','B']) * sxew_scale.loc[c,'Water Scale']
        intermediate11[intermediate11>sxew_cutoff_water] = sxew_cutoff_water
        op_sxew_water.loc[:,c] = (intermediate11 * \
         mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New SX-EW - Water
        intermediate12 = ((mine_eqn_coefficients.loc['SX-EW Water','A'] * mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Water','B']) * sxew_scale.loc[c,'Water Scale']).apply(\
         lambda x: x * new_sxew_dist.loc[:,c])
        intermediate12[intermediate12>sxew_cutoff_water] = sxew_cutoff_water
        new_sxew_water.loc[:,c] = (intermediate12 * \
         mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 

        country_impacts.loc[idx['20180101':,'Operating Traditional Mines'],idx[c,traci_impacts]] = op_mines_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating Traditional Mines'],idx[c,energy_impacts]] = op_mines_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating Traditional Mines'],idx[c,water_impacts]] = op_mines_water.loc[:,c].values
        country_impacts.loc[idx['20180101':,'Operating SX-EW Mines'],idx[c,traci_impacts]] = op_sxew_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating SX-EW Mines'],idx[c,energy_impacts]] = op_sxew_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating SX-EW Mines'],idx[c,water_impacts]] = op_sxew_water.loc[:,c].values

        country_impacts.loc[idx['20180101':,'New Traditional Mines'],idx[c,traci_impacts]] = new_mines_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'New Traditional Mines'],idx[c,energy_impacts]] = new_mines_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'New Traditional Mines'],idx[c,water_impacts]] = new_mines_water.loc[:,c].values
        country_impacts.loc[idx['20180101':,'New SX-EW Mines'],idx[c,traci_impacts]] = new_sxew_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'New SX-EW Mines'],idx[c,energy_impacts]] = new_sxew_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'New SX-EW Mines'],idx[c,water_impacts]] = new_sxew_water.loc[:,c].values

        country_quantities.loc[idx['20180101':, 'Operating Traditional Mines'],c] = mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].sum(axis=1).values
        country_quantities.loc[idx['20180101':, 'Operating SX-EW Mines'],c] = mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].sum(axis=1).values
        country_quantities.loc[idx['20180101':, 'New Traditional Mines'],c] = new_conc_dist.loc['20180101',c]*mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].sum(axis=1).values
        country_quantities.loc[idx['20180101':, 'New SX-EW Mines'],c] = new_sxew_dist.loc['20180101',c]*mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].sum(axis=1).values

        # Historical Mines
        country_quantities.loc[idx[:'20170101', 'Operating Traditional Mines'],c] = (sd_simulated.loc['19600101':'20170101','Concentrate production']\
            * conc_distribution.loc[c]).values
        country_quantities.loc[idx[:'20170101', 'Operating SX-EW Mines'],c] = (sd_simulated.loc['19600101':'20170101','SX-EW production']\
            * sxew_distribution.loc[c]).values
        country_impacts.loc[idx[:'20170101', 'Operating Traditional Mines'],idx[c,:]] = \
            (country_quantities.loc[idx[:'20170101', 'Operating Traditional Mines'],c].apply(lambda x: x*conc_impacts.loc[c,:])).values
        country_impacts.loc[idx[:'20170101', 'Operating SX-EW Mines'],idx[c,:]] = \
            (country_quantities.loc[idx[:'20170101', 'Operating SX-EW Mines'],c].apply(lambda x: x*sxew_impacts.loc[c,:])).values

        ###### Fabrication ######
        if c == 'China':
            country_impacts.loc[idx[:,'Fabrication'],idx[c,:]] = ((direct_melt_demand_cn.sum(axis=1))\
                .apply(lambda x: x * semis_impacts.loc['Metal working RoW (Conseq)', :])).values
            country_quantities.loc[idx[:,'Fabrication'],c] = (direct_melt_demand_cn.sum(axis=1)).values
        elif c == 'Europe' or c == 'North America':
            country_impacts.loc[idx[:,'Fabrication'],idx[c,:]] = ((direct_melt_demand_rw.sum(axis=1) * 
                mod_semis_distribution.loc[c].values)\
                .apply(lambda x: x * semis_impacts.loc['Metal working RER (Conseq)', :])).values
            country_quantities.loc[idx[:,'Fabrication'],c] = (direct_melt_demand_rw.sum(axis=1) * 
                mod_semis_distribution.loc[c].values).values
        else:
            country_impacts.loc[idx[:,'Fabrication'],idx[c,:]] = ((direct_melt_demand_rw.sum(axis=1) * 
                mod_semis_distribution.loc[c].values)\
                .apply(lambda x: x * semis_impacts.loc['Metal working RoW (Conseq)', :])).values
            country_quantities.loc[idx[:,'Fabrication'],c] = (direct_melt_demand_rw.sum(axis=1) * 
                mod_semis_distribution.loc[c].values).values

        ###### Direct Melt of Scrap and Other ######
        if c == 'China':
            country_impacts.loc[idx[:,'Direct Melt Scrap'], idx[c,:]] = (\
                direct_melt_demand_cn.loc[:, 'Yellow_Brass':'Red_Brass'].sum(axis=1).apply(lambda x : x * scrap_impacts.loc['Brass', :]) + \
                direct_melt_demand_cn.loc[:, 'Mn_Bronze':'Cartridge'].sum(axis=1).apply(lambda x : x * scrap_impacts.loc['Low Grade', :]) + \
                direct_melt_demand_cn.loc[:, 'No.1'].apply(lambda x : x * scrap_impacts.loc['No1', :]) + \
                direct_melt_demand_cn.loc[:, 'No.2'].apply(lambda x: x * scrap_impacts.loc['No2', :])).values
            country_impacts.loc[idx[:,'Direct Melt Other'], idx[c,:]] = ((\
                direct_melt_demand_cn.loc[:,:'Ref_Fe'].sum(axis=1)).apply(lambda x: x * scrap_impacts.loc['No1', :])).values # assuming that all the refined metals can be melted with the same impacts as No1 scrap
            country_quantities.loc[idx[:,'Direct Melt Scrap'],c] = (direct_melt_demand_cn.loc[:, 'Yellow_Brass':].sum(axis=1)).values
            country_quantities.loc[idx[:,'Direct Melt Other'],c] = (direct_melt_demand_cn.loc[:, :'Ref_Fe'].sum(axis=1)).values
        else:
            country_impacts.loc[idx[:,'Direct Melt Scrap'], idx[c,:]] = ((\
                direct_melt_demand_rw.loc[:, 'Yellow_Brass':'Red_Brass'].sum(axis=1).apply(lambda x : x * scrap_impacts.loc['Brass', :]) + \
                direct_melt_demand_rw.loc[:, 'Mn_Bronze':'Cartridge'].sum(axis=1).apply(lambda x : x * scrap_impacts.loc['Low Grade', :]) + \
                direct_melt_demand_rw.loc[:, 'No.1'].apply(lambda x : x * scrap_impacts.loc['No1', :]) + \
                direct_melt_demand_rw.loc[:, 'No.2'].apply(lambda x: x * scrap_impacts.loc['No2', :])) *\
                mod_semis_distribution.loc[c].values[0]).values
            country_impacts.loc[idx[:,'Direct Melt Other'], idx[c,:]] = ((\
                direct_melt_demand_rw.loc[:,:'Ref_Fe'].sum(axis=1)).apply(lambda x: x * scrap_impacts.loc['No1', :]) *
                mod_semis_distribution.loc[c].values[0]).values # assuming that all the refined metals can be melted with the same impacts as No1 scrap
            country_quantities.loc[idx[:,'Direct Melt Scrap'],c] = (direct_melt_demand_rw.loc[:, 'Yellow_Brass':].sum(axis=1)*mod_semis_distribution.loc[c].values[0]).values
            country_quantities.loc[idx[:,'Direct Melt Other'],c] = (direct_melt_demand_rw.loc[:, :'Ref_Fe'].sum(axis=1)*mod_semis_distribution.loc[c].values[0]).values

        ###### Refineries ###### 
        if c == 'China':
            country_impacts.loc[idx[:,'Primary Refinery'],idx[c,:]] = (sd_simulated_cn.loc['19600101':, 'Primary refining production']\
                .apply(lambda x: x * refinery_primary_impacts.loc[c, :])).values
            country_impacts.loc[idx[:,'Secondary Refinery'],idx[c,:]] = (sd_simulated_cn.loc['19600101':, 'Secondary refining production']\
                .apply(lambda x: x * refinery_secondary_impacts.loc[c, :])).values
            country_quantities.loc[idx[:,'Primary Refinery'],c] = (sd_simulated_cn.loc['19600101':, 'Primary refining production']).values
            country_quantities.loc[idx[:,'Secondary Refinery'],c] = (sd_simulated_cn.loc['19600101':, 'Secondary refining production']).values
        else:
            country_impacts.loc[idx[:,'Primary Refinery'],idx[c,:]] = ((sd_simulated_rw.loc['19600101':, 'Primary refining production'] * \
                mod_refinery_distribution.loc[c]).apply(lambda x: x * refinery_primary_impacts.loc[c, :])).values
            country_impacts.loc[idx[:,'Secondary Refinery'],idx[c,:]] = ((sd_simulated_rw.loc['19600101':, 'Secondary refining production'] * \
                mod_refinery_distribution.loc[c]).apply(lambda x: x * refinery_secondary_impacts.loc[c, :])).values
            country_quantities.loc[idx[:,'Primary Refinery'],c] = (sd_simulated_rw.loc['19600101':, 'Primary refining production']\
                * mod_refinery_distribution.loc[c]).values
            country_quantities.loc[idx[:,'Secondary Refinery'],c] = (sd_simulated_rw.loc['19600101':, 'Secondary refining production']\
                * mod_refinery_distribution.loc[c]).values

        ###### Refined metals' contributions ###### 
        if c == 'China':
            for m in refined_metal_impacts.index:
                country_impacts.loc[idx[:,'Refined Metals'],idx[c,:]] += ((direct_melt_demand_cn.loc[:,m])\
                    .apply(lambda x: x * refined_metal_impacts.loc[m, :])).values
            country_quantities.loc[idx[:,'Refined Metals'],c] = (\
                direct_melt_demand_cn.loc[:,refined_metal_impacts.index].sum(axis=1)).values
        else:
            for m in refined_metal_impacts.index:
                country_impacts.loc[idx[:,'Refined Metals'],idx[c,:]] += ((direct_melt_demand_rw.loc[:,m] * 
                    mod_semis_distribution.loc[c].values)\
                    .apply(lambda x: x * refined_metal_impacts.loc[m, :])).values
            country_quantities.loc[idx[:,'Refined Metals'],c] = (\
                direct_melt_demand_rw.loc[:,refined_metal_impacts.index].sum(axis=1)*
                mod_semis_distribution.loc[c].values).values

    colsss = pd.MultiIndex.from_product([countries, ['Direct Melt Scrap', 'Direct Melt Other', 
                 'Fabrication', 'New SX-EW Mines',
                 'New Traditional Mines', 'Operating SX-EW Mines', 'Operating Traditional Mines', 
                 'Primary Refinery', 'Secondary Refinery', 'Refined Metals']])
    idxss = pd.MultiIndex.from_product([direct_melt_demand_cn.index, impacts])
    new_country_impacts = pd.DataFrame(0, index = idxss, columns = colsss)
    new_country_unit_impacts = pd.DataFrame(0, index = idxss, columns = colsss)

    for y in country_impacts.index.levels[0]:
        for c in countries:
            new_country_impacts.loc[idx[y,:],idx[c,:]] = country_impacts.loc[idx[y,:],idx[c,:]].transpose().values
    
    # Getting country unit impacts
    for i in impacts:
        country_unit_impacts.loc[:, idx[:, i]] = (country_impacts.loc[:, idx[:, i]] / country_quantities.values).fillna(0)

    # Updating the columns to have the correct units (all the impacts imported were in kg X/kg Cu and we're multiplying by kt)
    impacts_units = dict()
    for i in impacts:
        impacts_units.update({i: i.replace('kg','kt').replace('CTU','kCTU').replace('MJ','GJ').replace('m3','ML')
                         .replace('Total','Total energy')})
#     print('Done')
    return [new_country_impacts.rename(index = impacts_units), 
            country_unit_impacts.rename(columns = impacts_units),
            country_quantities.rename(columns = impacts_units)]

def calculate_impacts_lin(sd_simulated, mine_life_stats_all_tp):
    '''Calculates environmental impacts based on production, returns total impacts, unit impacts, and quantities'''
    sd_simulated = sd_simulated.loc['19600101':].copy()
    idx = pd.IndexSlice
    mine_life_stats_all = mine_life_stats_all_tp.transpose()
    sxew_id_operating = [i for i in mine_life_stats_all.columns.levels[0] if 
                    round(mine_life_stats_all.loc[:, idx[i,'Recovered metal production (kt)']].mean()/
                     mine_life_stats_all.loc[:, idx[i,'Paid metal production (kt)']].mean(),6) <= 1
                     and 'Inc_' not in str(i)]
    conc_id_operating = [i for i in mine_life_stats_all.columns.levels[0] if 
                    i not in sxew_id_operating and 'Inc_' not in str(i)]
    sxew_id_new = [i for i in mine_life_stats_all.columns.levels[0] if 
                    round(mine_life_stats_all.loc[:, idx[i,'Recovered metal production (kt)']].mean()/
                     mine_life_stats_all.loc[:, idx[i,'Paid metal production (kt)']].mean(),6) <= 1
                     and 'Inc_' in str(i)]
    conc_id_new = [i for i in mine_life_stats_all.columns.levels[0] if 
                    i not in sxew_id_new and 'Inc_' in str(i)]
    id_new = [i for i in mine_life_stats_all.columns.levels[0] if 'Inc_' in str(i)]
    id_op = [i for i in mine_life_stats_all.columns.levels[0] if 'Inc_' not in str(i)]
    mine_life_stats_panel_operating = mine_life_stats_all.loc[:,idx[id_op,:]].copy()
    mine_life_stats_panel_new = mine_life_stats_all.loc[:,idx[id_new,:]].copy()
    
    # Loading data:
    import os
    directory = os.path.dirname(os.path.dirname(os.getcwd()))
    semis_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Semis Distribution', index_col=0).loc[:'Latin America',:'Fraction Semis Production']
    semis_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Semis', index_col=0)
    mod_semis_distribution = semis_distribution / (1 - semis_distribution.loc['China'])
    scrap_distribution = pd.Series([0.25,0.251969,0.348425,0.149606], index=['No1','No2','Brass','Low Grade'])
    
    refinery_primary_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Primary Refinery Ecoinvent', index_col=0).iloc[:,1:]
    refinery_secondary_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Scrap Ecoinvent', index_col=0, skiprows = 12).loc[:,:'Water use (m3)']
    refinery_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Primary Refinery Ecoinvent', index_col=0).iloc[:,0]
    mod_refinery_distribution = refinery_distribution / (1 - refinery_distribution.loc['China'])

    conc_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    conc_mult = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Mult', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    conc_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).iloc[:,0]
    conc_scale = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).loc[:,'TRACI Scale':]
    sxew_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Ecoinvent', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    sxew_mult = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Mult', index_col=0).loc[:,'Ozone depletion (kg CFC-kk)':'Water use (m3)']
    sxew_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Ecoinvent', index_col=0).iloc[:,0]/100
    sxew_scale = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'SX-EW Ecoinvent', index_col=0).loc[:,'TRACI Scale':]
    mining_evolution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Conc Ecoinvent', index_col=0).loc[:,['Fraction','Identified Resources','Undiscovered Resources']]
    mine_eqn_coefficients = pd.DataFrame([[1.5777, -0.626], [2.0621, -1.208], [15.697, -0.573], [36.529, -0.351], 
                                          [0.073378, -0.094], [0.041232, -0.34]], 
                                         index = ['Conc TRACI', 'SX-EW TRACI','Conc Energy','SX-EW Energy',
                                                  'Conc Water','SX-EW Water'], columns = ['A','B'])
    mine_eqn_coefficients.index.name = '(A*ore_grade^B) * scale'

    scrap_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Scrap Ecoinvent', index_col=0, skiprows = 1).loc[:'No1',:'Water use (m3)']
    # direct_melt_distribution = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Semis Distribution', index_col=0).loc[:'Latin America',:'Fraction Semis Production']
    refined_metal_impacts = pd.read_excel(directory+'/Data/SimaPro Data.xlsx', sheet_name = 'Ref_Metals', index_col=0).iloc[:,2:]
    refined_metal_impacts.index = ['Ref_Zn','Ref_Pb','Ref_Sn','Ref_Ni','Ref_Al','Ref_Mn','Ref_Fe']
    mine_regionality = pd.read_excel(directory+'/Data/primary supply/Operating mine pool - countries.xlsx', index_col=0)

    # Mining evolution
    new_conc_dist = pd.DataFrame(np.array((list(mining_evolution.loc[:,'Fraction'])*23)).reshape(23, 7), index = mine_life_stats_panel_operating.index, columns = mining_evolution.index)
    new_sxew_dist = pd.DataFrame(np.array((list(sxew_distribution)*23)).reshape(23, 7), index = mine_life_stats_panel_operating.index, columns = sxew_distribution.index)

    # Other initialization
    # Semis first (excluding refined metals, those will be last so we have that year's copper impacts)
    countries = refinery_primary_impacts.index
    impacts = refinery_primary_impacts.columns
    traci_impacts = impacts[:10]
    energy_impacts = impacts[10:-1]
    water_impacts = list([impacts[-1]])
    colss = pd.MultiIndex.from_product([countries, impacts])
    idxs = pd.MultiIndex.from_product([sd_simulated.loc['19600101':].index, \
                ['Direct Melt Scrap', 'Direct Melt Other', 'Fabrication', 'New SX-EW Mines',
                 'New Traditional Mines', 'Operating SX-EW Mines', 'Operating Traditional Mines', 
                 'Primary Refinery', 'Secondary Refinery', 'Refined Metals']])
    country_impacts = pd.DataFrame(0, index = idxs, columns = colss)
    country_unit_impacts = pd.DataFrame(0, index = idxs, columns = colss)
    country_quantities = pd.DataFrame(0, index = idxs, columns = countries)

    # Mine and SX-EW preliminary calculations
    op_mines_traci = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    new_mines_traci = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    op_sxew_traci = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    new_sxew_traci = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    op_mines_energy = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    new_mines_energy = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    op_sxew_energy = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    new_sxew_energy = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    op_mines_water = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    new_mines_water = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    op_sxew_water = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)
    new_sxew_water = pd.DataFrame(0, index = sd_simulated.loc['20180101':].index, columns = countries)

    direct_melt_demand_rw = sd_simulated.loc['19600101':,['Refined demand','Direct melt scrap']].copy()
    sd_simulated_rw = sd_simulated.loc['19600101':].copy()
        
    for c in countries:
        ####### Mining #######
        conc_op_region_id = [i for i in conc_id_operating if c in mine_regionality.loc[i,'Region']]
        sxew_op_region_id = [i for i in sxew_id_operating if c in mine_regionality.loc[i,'Region']]
        conc_cutoff_traci = 10
        sxew_cutoff_traci = 28
        conc_cutoff_energy = 150
        sxew_cutoff_energy = 450
        conc_cutoff_water = 0.27
        sxew_cutoff_water = 0.1

        # Operating Concentrate Mines - TRACI
        intermediate1 = (mine_eqn_coefficients.loc['Conc TRACI','A'] * mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc TRACI','B']) * conc_scale.loc[c,'TRACI Scale']
        intermediate1[intermediate1>conc_cutoff_traci] = conc_cutoff_traci
        op_mines_traci.loc[:,c] = (intermediate1 * \
         mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New Concentrate Mines - TRACI
        intermediate2 = ((mine_eqn_coefficients.loc['Conc TRACI','A'] * mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc TRACI','B']) * conc_scale.loc[c,'TRACI Scale']).apply(\
         lambda x: x * new_conc_dist.loc[:,c])
        intermediate2[intermediate2>conc_cutoff_traci] = conc_cutoff_traci
        new_mines_traci.loc[:,c] = (intermediate2 * \
         mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating SX-EW - TRACI
        intermediate3 = (mine_eqn_coefficients.loc['SX-EW TRACI','A'] * mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW TRACI','B']) * sxew_scale.loc[c,'TRACI Scale']
        intermediate3[intermediate3>sxew_cutoff_traci] = sxew_cutoff_traci
        op_sxew_traci.loc[:,c] = (intermediate3 * \
         mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New SX-EW - TRACI
        intermediate4 = ((mine_eqn_coefficients.loc['SX-EW TRACI','A'] * mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW TRACI','B']) * sxew_scale.loc[c,'TRACI Scale']).apply(\
         lambda x: x * new_sxew_dist.loc[:,c])
        intermediate4[intermediate4>sxew_cutoff_traci] = sxew_cutoff_traci
        new_sxew_traci.loc[:,c] = (intermediate4 * \
         mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating Concentrate Mines - Energy
        intermediate5 = (mine_eqn_coefficients.loc['Conc Energy','A'] * mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Energy','B'] ) * conc_scale.loc[c,'Energy Scale']
        intermediate5[intermediate5>conc_cutoff_energy] = conc_cutoff_energy
        op_mines_energy.loc[:,c] = (intermediate5 * \
         mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New Concentrate Mines - Energy
        intermediate6 = ((mine_eqn_coefficients.loc['Conc Energy','A'] * mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Energy','B']) * conc_scale.loc[c,'Energy Scale']).apply(\
         lambda x: x * new_conc_dist.loc[:,c])
        intermediate6[intermediate6>conc_cutoff_energy] = conc_cutoff_energy
        new_mines_energy.loc[:,c] = (intermediate6 * \
         mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating SX-EW - Energy
        intermediate7 = (mine_eqn_coefficients.loc['SX-EW Energy','A'] * mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Energy','B']) * sxew_scale.loc[c,'Energy Scale']
        intermediate7[intermediate7>sxew_cutoff_energy] = sxew_cutoff_energy
        op_sxew_energy.loc[:,c] = (intermediate7 * \
         mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New SX-EW - Energy
        intermediate8 = ((mine_eqn_coefficients.loc['SX-EW Energy','A'] * mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Energy','B']) * sxew_scale.loc[c,'Energy Scale']).apply(\
         lambda x: x * new_sxew_dist.loc[:,c])
        intermediate8[intermediate8>sxew_cutoff_energy] = sxew_cutoff_energy
        new_sxew_energy.loc[:,c] = (intermediate8 * \
         mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating Concentrate Mines - Water
        intermediate9 = (mine_eqn_coefficients.loc['Conc Water','A'] * mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Water','B'] ) * conc_scale.loc[c,'Water Scale']
        intermediate9[intermediate9>conc_cutoff_water] = conc_cutoff_water
        op_mines_water.loc[:,c] = (intermediate9 * \
         mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New Concentrate Mines - Water
        intermediate10 = ((mine_eqn_coefficients.loc['Conc Water','A'] * mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['Conc Water','B']) * conc_scale.loc[c,'Water Scale']).apply(\
         lambda x: x * new_conc_dist.loc[:,c])
        intermediate10[intermediate10>conc_cutoff_water] = conc_cutoff_water
        new_mines_water.loc[:,c] = (intermediate10 * \
         mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 
        # Operating SX-EW - Water
        intermediate11 = (mine_eqn_coefficients.loc['SX-EW Water','A'] * mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Water','B']) * sxew_scale.loc[c,'Water Scale']
        intermediate11[intermediate11>sxew_cutoff_water] = sxew_cutoff_water
        op_sxew_water.loc[:,c] = (intermediate11 * \
         mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].values).sum(axis=1).values
        # New SX-EW - Water
        intermediate12 = ((mine_eqn_coefficients.loc['SX-EW Water','A'] * mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Head grade (%)']]**\
         mine_eqn_coefficients.loc['SX-EW Water','B']) * sxew_scale.loc[c,'Water Scale']).apply(\
         lambda x: x * new_sxew_dist.loc[:,c])
        intermediate12[intermediate12>sxew_cutoff_water] = sxew_cutoff_water
        new_sxew_water.loc[:,c] = (intermediate12 * \
         mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].values).sum(axis=1).values 

        country_impacts.loc[idx['20180101':,'Operating Traditional Mines'],idx[c,traci_impacts]] = op_mines_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating Traditional Mines'],idx[c,energy_impacts]] = op_mines_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating Traditional Mines'],idx[c,water_impacts]] = op_mines_water.loc[:,c].values
        country_impacts.loc[idx['20180101':,'Operating SX-EW Mines'],idx[c,traci_impacts]] = op_sxew_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating SX-EW Mines'],idx[c,energy_impacts]] = op_sxew_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'Operating SX-EW Mines'],idx[c,water_impacts]] = op_sxew_water.loc[:,c].values

        country_impacts.loc[idx['20180101':,'New Traditional Mines'],idx[c,traci_impacts]] = new_mines_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'New Traditional Mines'],idx[c,energy_impacts]] = new_mines_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'New Traditional Mines'],idx[c,water_impacts]] = new_mines_water.loc[:,c].values
        country_impacts.loc[idx['20180101':,'New SX-EW Mines'],idx[c,traci_impacts]] = new_sxew_traci.loc[:,c].apply(lambda x: x*conc_mult.loc[c,traci_impacts]).values
        country_impacts.loc[idx['20180101':,'New SX-EW Mines'],idx[c,energy_impacts]] = new_sxew_energy.loc[:,c].apply(lambda x: x*conc_mult.loc[c,energy_impacts]).values
        country_impacts.loc[idx['20180101':,'New SX-EW Mines'],idx[c,water_impacts]] = new_sxew_water.loc[:,c].values

        country_quantities.loc[idx['20180101':, 'Operating Traditional Mines'],c] = mine_life_stats_panel_operating.loc[:,idx[conc_op_region_id,'Recovered metal production (kt)']].sum(axis=1).values
        country_quantities.loc[idx['20180101':, 'Operating SX-EW Mines'],c] = mine_life_stats_panel_operating.loc[:,idx[sxew_op_region_id,'Recovered metal production (kt)']].sum(axis=1).values
        country_quantities.loc[idx['20180101':, 'New Traditional Mines'],c] = mine_life_stats_panel_new.loc[:,idx[conc_id_new,'Recovered metal production (kt)']].sum(axis=1).values
        country_quantities.loc[idx['20180101':, 'New SX-EW Mines'],c] = mine_life_stats_panel_new.loc[:,idx[sxew_id_new,'Recovered metal production (kt)']].sum(axis=1).values

        # Historical Mines
        country_quantities.loc[idx[:'20170101', 'Operating Traditional Mines'],c] = (sd_simulated.loc['19600101':'20170101','Concentrate production']\
            * conc_distribution.loc[c]).values
        country_quantities.loc[idx[:'20170101', 'Operating SX-EW Mines'],c] = (sd_simulated.loc['19600101':'20170101','SX-EW production']\
            * sxew_distribution.loc[c]).values
        country_impacts.loc[idx[:'20170101', 'Operating Traditional Mines'],idx[c,:]] = \
            (country_quantities.loc[idx[:'20170101', 'Operating Traditional Mines'],c].apply(lambda x: x*conc_impacts.loc[c,:])).values
        country_impacts.loc[idx[:'20170101', 'Operating SX-EW Mines'],idx[c,:]] = \
            (country_quantities.loc[idx[:'20170101', 'Operating SX-EW Mines'],c].apply(lambda x: x*sxew_impacts.loc[c,:])).values

        ###### Fabrication ######
        if c == 'Europe' or c == 'North America':
            country_impacts.loc[idx[:,'Fabrication'],idx[c,:]] = ((direct_melt_demand_rw.sum(axis=1) * 
                semis_distribution.loc[c].values)\
                .apply(lambda x: x * semis_impacts.loc['Metal working RER (Conseq)', :])).values
            country_quantities.loc[idx[:,'Fabrication'],c] = (direct_melt_demand_rw.sum(axis=1) * 
                semis_distribution.loc[c].values).values
        else:
            country_impacts.loc[idx[:,'Fabrication'],idx[c,:]] = ((direct_melt_demand_rw.sum(axis=1) * 
                semis_distribution.loc[c].values)\
                .apply(lambda x: x * semis_impacts.loc['Metal working RoW (Conseq)', :])).values
            country_quantities.loc[idx[:,'Fabrication'],c] = (direct_melt_demand_rw.sum(axis=1) * 
                semis_distribution.loc[c].values).values

        ###### Direct Melt of Scrap and Other ######
        country_impacts.loc[idx[:,'Direct Melt Scrap'], idx[c,:]] = ((\
            (scrap_distribution.loc['Brass']*sd_simulated.loc[:,'Direct melt scrap']).apply(lambda x : x * scrap_impacts.loc['Brass', :]) + \
            (scrap_distribution.loc['Low Grade']*sd_simulated.loc[:,'Direct melt scrap']).apply(lambda x : x * scrap_impacts.loc['Low Grade', :]) + \
            (scrap_distribution.loc['No1']*sd_simulated.loc[:,'Direct melt scrap']).apply(lambda x : x * scrap_impacts.loc['No1', :]) + \
            (scrap_distribution.loc['No2']*sd_simulated.loc[:,'Direct melt scrap']).apply(lambda x: x * scrap_impacts.loc['No2', :])) *\
            semis_distribution.loc[c].values[0]).values
        country_impacts.loc[idx[:,'Direct Melt Other'], idx[c,:]] = ((\
            sd_simulated.loc[:,'Refined demand']).apply(lambda x: x * scrap_impacts.loc['No1', :]) *
            semis_distribution.loc[c].values[0]).values # assuming that all the refined metals can be melted with the same impacts as No1 scrap
        country_quantities.loc[idx[:,'Direct Melt Scrap'],c] = (sd_simulated.loc[:,'Direct melt scrap']*semis_distribution.loc[c].values[0]).values
        country_quantities.loc[idx[:,'Direct Melt Other'],c] = (sd_simulated.loc[:,'Refined demand']*semis_distribution.loc[c].values[0]).values

        ###### Refineries ###### 
        country_impacts.loc[idx[:,'Primary Refinery'],idx[c,:]] = ((sd_simulated_rw.loc['19600101':, 'Primary refining production'] * \
            refinery_distribution.loc[c]).apply(lambda x: x * refinery_primary_impacts.loc[c, :])).values
        country_impacts.loc[idx[:,'Secondary Refinery'],idx[c,:]] = ((sd_simulated_rw.loc['19600101':, 'Secondary refining production'] * \
            refinery_distribution.loc[c]).apply(lambda x: x * refinery_secondary_impacts.loc[c, :])).values
        country_quantities.loc[idx[:,'Primary Refinery'],c] = (sd_simulated_rw.loc['19600101':, 'Primary refining production']\
            * refinery_distribution.loc[c]).values
        country_quantities.loc[idx[:,'Secondary Refinery'],c] = (sd_simulated_rw.loc['19600101':, 'Secondary refining production']\
            * refinery_distribution.loc[c]).values

        ###### Refined metals' contributions ###### 
#         if c == 'China':
#             for m in refined_metal_impacts.index:
#                 country_impacts.loc[idx[:,'Refined Metals'],idx[c,:]] += ((direct_melt_demand_cn.loc[:,m])\
#                     .apply(lambda x: x * refined_metal_impacts.loc[m, :])).values
#             country_quantities.loc[idx[:,'Refined Metals'],c] = (\
#                 direct_melt_demand_cn.loc[:,refined_metal_impacts.index].sum(axis=1)).values
#         else:
#             for m in refined_metal_impacts.index:
#                 country_impacts.loc[idx[:,'Refined Metals'],idx[c,:]] += ((direct_melt_demand_rw.loc[:,m] * 
#                     mod_semis_distribution.loc[c].values)\
#                     .apply(lambda x: x * refined_metal_impacts.loc[m, :])).values
#             country_quantities.loc[idx[:,'Refined Metals'],c] = (\
#                 direct_melt_demand_rw.loc[:,refined_metal_impacts.index].sum(axis=1)*
#                 mod_semis_distribution.loc[c].values).values

    colsss = pd.MultiIndex.from_product([countries, ['Direct Melt Scrap', 'Direct Melt Other', 
                 'Fabrication', 'New SX-EW Mines',
                 'New Traditional Mines', 'Operating SX-EW Mines', 'Operating Traditional Mines', 
                 'Primary Refinery', 'Secondary Refinery', 'Refined Metals']])
    idxss = pd.MultiIndex.from_product([sd_simulated.loc['19600101':].index, impacts])
    new_country_impacts = pd.DataFrame(0, index = idxss, columns = colsss)
    new_country_unit_impacts = pd.DataFrame(0, index = idxss, columns = colsss)

    for y in country_impacts.index.levels[0]:
        for c in countries:
            new_country_impacts.loc[idx[y,:],idx[c,:]] = country_impacts.loc[idx[y,:],idx[c,:]].transpose().values
    
    # Getting country unit impacts
    for i in impacts:
        country_unit_impacts.loc[:, idx[:, i]] = (country_impacts.loc[:, idx[:, i]] / country_quantities.values).fillna(0)

    # Updating the columns to have the correct units (all the impacts imported were in kg X/kg Cu and we're multiplying by kt)
    impacts_units = dict()
    for i in impacts:
        impacts_units.update({i: i.replace('kg','kt').replace('CTU','kCTU').replace('MJ','GJ').replace('m3','ML')
                         .replace('Total','Total energy')})
#     print('Done')
    return [new_country_impacts.rename(index = impacts_units), 
            country_unit_impacts.rename(columns = impacts_units),
            country_quantities.rename(columns = impacts_units)]

def pct_inc(impact_change,metric,baseline):
    '''Returns, global, China, and RoW responses'''
    ph_series1 = impact_change.groupby(axis=1,level=1).sum().loc[idx[:,metric],:].reset_index(level=1,drop=True)
    ph_series2 = impact_change.loc[idx[:,metric],idx['China', :]].reset_index(level=1,drop=True).T.reset_index(level=0,drop=True).T
    ph_series3 = impact_change.loc[idx[:,metric],idx[['Oceania','Africa','Europe','North America','Latin America','Other Asia'],:]].groupby(axis=1,level=1).sum().reset_index(level=1,drop=True)
    ph_series1 = ph_series1.cumsum().loc[:,categories]/total_impacts_pir_baseline.loc[idx['20200101',metric]].sum()*100
    ph_series2 = ph_series2.cumsum().loc[:,categories]/total_impacts_pir_baseline.loc[idx['20200101',metric]].sum()*100
    ph_series3 = ph_series3.cumsum().loc[:,categories]/total_impacts_pir_baseline.loc[idx['20200101',metric]].sum()*100
    return ph_series1,ph_series2,ph_series3

def abs_inc(impact_change,metric,baseline):
    '''Returns, global, China, and RoW responses'''
    ph_series1 = impact_change.groupby(axis=1,level=1).sum().loc[idx[:,metric],:].reset_index(level=1,drop=True)
    ph_series2 = impact_change.loc[idx[:,metric],idx['China', :]].reset_index(level=1,drop=True).T.reset_index(level=0,drop=True).T
    ph_series3 = impact_change.loc[idx[:,metric],idx[['Oceania','Africa','Europe','North America','Latin America','Other Asia'],:]].groupby(axis=1,level=1).sum().reset_index(level=1,drop=True)
    ph_series1 = ph_series1.cumsum().loc[:,categories]
    ph_series2 = ph_series2.cumsum().loc[:,categories]
    ph_series3 = ph_series3.cumsum().loc[:,categories]
    return ph_series1,ph_series2,ph_series3

def pct_inc_lin(impact_change,metric,baseline):
    '''Returns global response as fraction of 2020'''
    ph_series1 = impact_change.groupby(axis=1,level=1).sum().loc[idx[:,metric],:].reset_index(level=1,drop=True)
    ph_series1 = ph_series1.cumsum().loc[:,categories]/baseline.loc[idx['20200101',metric]].sum()*100
    return ph_series1

def abs_inc_lin(impact_change,metric,baseline):
    '''Returns absolute global response'''
    ph_series1 = impact_change.groupby(axis=1,level=1).sum().loc[idx[:,metric],:].reset_index(level=1,drop=True)
    ph_series1 = ph_series1.cumsum().loc[:,categories]
    return ph_series1

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)
