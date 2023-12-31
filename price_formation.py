import pandas as pd
import numpy as np


def cathode_price_predict(cathode_price_last, cathode_balance, cathode_sd_elas):
    cathode_price_next=cathode_price_last+cathode_sd_elas*(cathode_balance)
    if cathode_price_next < 2000:
        cathode_price_next=2000
    elif cathode_price_next > 10000:
        cathode_price_next=10000
    return cathode_price_next



def tcrc_predict(tcrc_last, conc_balance, conc_sd_elas):
    tcrc_next=tcrc_last+conc_sd_elas*(conc_balance)
    if tcrc_next < 280:
        tcrc_next=280
    elif tcrc_next > 3000: #original: 1000, attempting 3000 for 100% shocks
        tcrc_next=3000
    return tcrc_next


def sp2_predict(sp2_last, scrap_balance, scrap_sd_elas, cathode_diff, 
                cathode_sp2_elas):
    sp2_next=sp2_last+scrap_sd_elas*scrap_balance+cathode_sp2_elas*cathode_diff
    if sp2_next < 300:
        sp2_next=300
    elif sp2_next > 5000: # original: 1250, 5000 seems to work for 100% shocks
        sp2_next=5000
    return sp2_next


