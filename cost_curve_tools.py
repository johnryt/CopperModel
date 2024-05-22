import pandas as pd
import numpy as np
idx = pd.IndexSlice


def cost_curve_df(cc_year, mine_life_stats_panel_all):
    cost_components=['Minsite cost ($/tonne paid metal)', 'TCRC ($/tonne paid metal)', 
                     'Frieght ($/tonne paid metal)', 'Royalty ($/tonne paid metal)']

    cc_t=datetime(cc_year, 1, 1)
    cost_components_of_year=mine_life_stats_panel_all.loc[cc_t, idx[:, cost_components]]
    cost_total_of_year=cost_components_of_year.groupby(level=0).sum()
    cost_orders=cost_total_of_year[cost_total_of_year>0].sort_values(ascending=True).index

    mine_index=mine_life_stats_panel_all.columns.get_level_values(0).unique()
    metal_prod=pd.Series(mine_life_stats_panel_all.loc[cc_t, idx[:, 'Recovered metal production (kt)']].values, 
                         index=mine_index).loc[cost_orders]

    cc_df=pd.DataFrame({'Total cost': cost_total_of_year.loc[cost_orders], 'Cumulative production': metal_prod.cumsum()})
    cc_df=pd.concat([cc_df, cost_components_of_year.unstack().loc[cost_orders,:]], axis=1)
    
    # Remove the 2018 Very high cost 'Outlier'
    if cc_year==2018:
        cc_df=cc_df.iloc[:-13, :]
    
    return cc_df


def cost_percentile(cc_df, percentile):
    prod_at_percentile=cc_df.loc[:, 'Cumulative production'].iloc[-1]*percentile/100

    cost_at_percentile=cc_df.iloc[(cc_df['Cumulative production']-prod_at_percentile)\
                                  .abs().argsort()[:2]].loc[:, 'Total cost'].mean()
    return cost_at_percentile


def cost_percentile_series(mine_life_stats_panel_all, percentile):
    cost_series=pd.Series(0, index=np.arange(2018, 2041))
    for year_i in np.arange(2018, 2041):
        cc_df_year_i=cost_curve_df(year_i, mine_life_stats_panel_all)
        cost_at_percentile_year_i=cost_percentile(cc_df_year_i, percentile)
        cost_series.loc[year_i]=cost_at_percentile_year_i
    return cost_series