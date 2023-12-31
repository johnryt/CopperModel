import pandas as pd


def cu_growth_cal(TCRC_growth, TCRC_elas):
    growth=TCRC_growth**(TCRC_elas)
    return growth


def sec_ratio_growth_cal(TCRC_growth, spread_growth, TCRC_elas, spread_elas):
    growth=TCRC_growth**(TCRC_elas)*spread_growth**(spread_elas)
    return growth


def pri_cap_growth_cal(concentrate_growth):
    growth=concentrate_growth
    return growth


def sec_cap_growth_cal(mine_growth, scrap_growth, sec_coef=0):
    growth=mine_growth*(1-sec_coef)+scrap_growth*sec_coef
    return growth


def simulate_refinery_production_oneyear(year_i, tcrc_series, sp2_series, 
                                         pri_cap_growth_series, sec_cap_growth_series,
                                         ref_stats, ref_hyper_param, 
                                         sec_coef=0, growth_lag=1):
    
    ref_stats_next=pd.Series(0, index=ref_stats.columns)
    
    pri_CU_TCRC_elas=ref_hyper_param.loc['pri CU TCRC elas', 'Value']
    sec_CU_TCRC_elas=ref_hyper_param.loc['sec CU TCRC elas', 'Value']
    sec_ratio_TCRC_elas=ref_hyper_param.loc['sec ratio TCRC elas', 'Value']
    sec_ratio_SP2_elas=ref_hyper_param.loc['sec ratio SP2 elas', 'Value']
    
    t=pd.datetime(year_i, 1, 1)
    t_lag_1=pd.datetime(year_i-1, 1, 1)
    t_lag_2=pd.datetime(year_i-2, 1, 1)
    
    tcrc_growth=tcrc_series.loc[t]/tcrc_series.loc[t_lag_1]
    sp2_growth=sp2_series.loc[t]/sp2_series.loc[t_lag_1]
    if growth_lag == 1:
        pri_growth=pri_cap_growth_series.loc[t_lag_1]/pri_cap_growth_series.loc[t_lag_2]
        sec_growth=sec_cap_growth_series.loc[year_i-1]/sec_cap_growth_series.loc[year_i-2]
    else:
        pri_growth=pri_cap_growth_series.loc[t]/pri_cap_growth_series.loc[t_lag_1]
        sec_growth=sec_cap_growth_series.loc[year_i]/sec_cap_growth_series.loc[year_i-1]
    
    pri_cu_growth=cu_growth_cal(tcrc_growth, pri_CU_TCRC_elas)
    sec_cu_growth=cu_growth_cal(tcrc_growth, sec_CU_TCRC_elas)
    ratio_growth=sec_ratio_growth_cal(tcrc_growth, sp2_growth, sec_ratio_TCRC_elas, sec_ratio_SP2_elas)
    pri_cap_growth=pri_cap_growth_cal(pri_growth)
    sec_cap_growth=sec_cap_growth_cal(pri_growth, sec_growth, sec_coef)

    ref_stats_next.loc['Primary CU']=ref_stats.loc[t_lag_1, 'Primary CU']*pri_cu_growth
    ref_stats_next.loc['Secondary CU']=ref_stats.loc[t_lag_1, 'Secondary CU']*sec_cu_growth
    ref_stats_next.loc['Secondary ratio']=ref_stats.loc[t_lag_1, 'Secondary ratio']*ratio_growth
    ref_stats_next.loc['Primary capacity']=ref_stats.loc[t_lag_1, 'Primary capacity']*pri_cap_growth
    ref_stats_next.loc['Secondary capacity']=ref_stats.loc[t_lag_1, 'Secondary capacity']*sec_cap_growth
    ref_stats_next.loc['Primary production']=ref_stats_next.loc['Primary capacity']*ref_stats_next.loc['Primary CU']\
    +ref_stats_next.loc['Secondary capacity']*ref_stats_next.loc['Secondary CU']*(1-ref_stats_next.loc['Secondary ratio'])
    ref_stats_next.loc['Secondary production']=\
    ref_stats_next.loc['Secondary capacity']*ref_stats_next.loc['Secondary CU']*ref_stats_next.loc['Secondary ratio']
    
    return ref_stats_next


def ref_stats_init(simulation_time, ref_hyper_param):
    ref_stats=pd.DataFrame(0, index=simulation_time, 
                           columns=['Primary capacity', 'Primary CU', 'Secondary capacity', 'Secondary CU',
                                    'Secondary ratio', 'Primary production', 'Secondary production'])
    pri_cap=ref_hyper_param.loc['pri cap', 'Value']
    pri_CU=ref_hyper_param.loc['pri CU', 'Value']
    sec_cap=ref_hyper_param.loc['sec cap', 'Value']
    sec_CU=ref_hyper_param.loc['sec CU', 'Value']
    sec_ratio=ref_hyper_param.loc['sec ratio', 'Value']
    
    ref_stats.loc['20180101', 'Primary capacity']=pri_cap
    ref_stats.loc['20180101', 'Primary CU']=pri_CU
    ref_stats.loc['20180101', 'Secondary capacity']=sec_cap
    ref_stats.loc['20180101', 'Secondary CU']=sec_CU
    ref_stats.loc['20180101', 'Secondary ratio']=sec_ratio
    
    ref_stats.loc['20180101', 'Primary production']=pri_cap*pri_CU+sec_cap*sec_CU*(1-sec_ratio)
    ref_stats.loc['20180101', 'Secondary production']=sec_cap*sec_CU*sec_ratio
    
    return ref_stats

