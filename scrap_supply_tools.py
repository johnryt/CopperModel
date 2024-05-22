import pandas as pd
import numpy as np
from scipy.stats import lognorm
from scipy.optimize import fsolve


# Lognormal mean function
def lognorm_mean(mu, sigma_to_mu, mean):
    sigma=sigma_to_mu*mu
    return np.exp(mu+sigma**2/2)-mean


# Calculate mu based on lognormal mean
def calc_mu(sigma_to_mu, mean, x0, **kwargs):
    return fsolve(lognorm_mean, x0=x0, args=(sigma_to_mu, mean), **kwargs).item()


# Build lifetime dataframe
def lifetime_df(product_lifetime, sigma_to_mu=0.1, mu_0=3):
    product_lifetime_df=pd.DataFrame(0, index=product_lifetime.index, columns=['Lifetime', 'mu', 'sigma'])
    product_lifetime_df.loc[:, 'Lifetime']=product_lifetime.values

    for p in product_lifetime.index:
        lifetime_mean=product_lifetime.loc[p]
        mu=calc_mu(sigma_to_mu, mean=lifetime_mean, x0=mu_0)
        product_lifetime_df.loc[p, 'mu']=mu
        product_lifetime_df.loc[p, 'sigma']=mu*sigma_to_mu
        
    return product_lifetime_df


# Return lifetime lognormal distribution of each product, in discrete frequencies
def lifetime_freq_df(product_lifetime_df):
    lifetime_freq = pd.DataFrame(0, index=product_lifetime_df.index, columns=np.arange(0,200,1))
    # First year reaching end of life is year 0, which 
    for p in lifetime_freq.index:
        mu = product_lifetime_df.loc[p, 'mu']
        sigma = product_lifetime_df.loc[p, 'sigma']
        freqs = np.diff(lognorm.cdf(np.arange(0,201,1), s=sigma, loc=0, scale=np.exp(mu)))
        lifetime_freq.loc[p, :] = freqs
    
    return lifetime_freq


# Simulate product reaching end of life, for all history
def product_reach_eol(use_product, product_lifetime_freq_df):
    product_eol = pd.DataFrame(0, index=use_product.index, columns=use_product.columns)
    for y in product_eol.index:
        for p in product_eol.columns:
            enter_use = use_product.loc[y, p]
            eol_freq = product_lifetime_freq_df.loc[p,:]
            eol_gen = enter_use * eol_freq
            eol_gen.index=eol_gen.index+y

            product_eol.loc[y:2018, p] += eol_gen.loc[y:2018]
            
    return product_eol


# Simulate product reaching end of life, for one year
def product_reach_eol_oneyear(year_i, use_product, product_lifetime_freq_df):
    first_use_year=use_product.index[0]
    product_eol=pd.Series(0, index=use_product.columns)
    for p in product_eol.index:
        for year_enter_use in np.arange(first_use_year, year_i+1):
            enter_use = use_product.loc[year_enter_use, p]
            eol_freq = product_lifetime_freq_df.loc[p,:]
            eol_gen = enter_use*eol_freq
            eol_gen.index=eol_gen.index+year_enter_use

            product_eol.loc[p] += eol_gen.loc[year_i]
            
    return product_eol


# Calculate waste collected by type, based on product reaching end of life
def waste_collected_oneyear(product_eol_year_i, product_to_waste_collectable, sort_eff, collect_rate):
    waste_eol=pd.Series(np.matmul(product_eol_year_i, product_to_waste_collectable), 
                        index=product_to_waste_collectable.columns)
    waste_collected=waste_eol.mul(sort_eff).mul(collect_rate)
    return waste_collected


# New scrap simulation of one year
def simulate_new_scrap_one_year(year_i, use_product, new_scrap_gen_ratio, 
                                product_to_waste_no_loss, sort_eff,
                                home_scrap_ratio=0.45, exchange_scrap_ratio=0.45):
    use_product_year_i=use_product.loc[year_i, :]
    new_scrap_product=use_product_year_i.mul(new_scrap_gen_ratio)
    exchange_scrap_product=new_scrap_product.mul(exchange_scrap_ratio)
    external_scrap_product=new_scrap_product.mul(1-home_scrap_ratio-exchange_scrap_ratio)
    external_scrap_waste=\
    pd.Series(np.matmul(external_scrap_product, product_to_waste_no_loss), 
              index=product_to_waste_no_loss.columns).mul(sort_eff)
    exchange_scrap_waste=\
    pd.Series(np.matmul(exchange_scrap_product, product_to_waste_no_loss), 
              index=product_to_waste_no_loss.columns)
    new_scrap_waste=external_scrap_waste+exchange_scrap_waste
    return new_scrap_waste