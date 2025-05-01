import pickle

import matplotlib.pylab as plt
import numpy as np
import pandas as pd

##CUTOFF

table = pd.read_csv('/Users/cretignier/Documents/THE/done_the_isochrones.csv',index_col=0)
table['dist'] = 1000/table['parallax']
#airmass = pd.read_pickle('/Users/cretignier/Documents/THE/Airmass_gr8.p')
table['index'] = np.arange(len(table))

def func_cutoff(table,cutoff,plot=True):
    plt.figure(figsize=(16,10))
    table2 = table.copy()
    count=0
    for kw in cutoff.keys():
        count+=1
        value = cutoff[kw]
        if kw[-1]=='<':
            mask = table2[kw[0:-1]]<value
        else:
            mask = table2[kw[0:-1]]>value
        plt.subplot(4,4,count)
        plt.title(kw+str(value))
        plt.hist(table2[kw[0:-1]],cumulative=True,bins=100)
        plt.axvline(x=value,label='%.0f / %.0f'%(sum(mask),len(mask)),color='k')
        if kw[-1]=='<':
            plt.axvspan(xmin=value,xmax=np.nanmax(table2[kw[0:-1]]),alpha=0.2,color='k')
        else:
            plt.axvspan(xmax=value,xmin=np.nanmin(table2[kw[0:-1]]),alpha=0.2,color='k')
        plt.legend()
        plt.grid()
        plt.xlim(np.nanmin(table2[kw[0:-1]]),np.nanmax(table2[kw[0:-1]]))
        table2 = table2[mask].reset_index(drop=True)
    plt.subplots_adjust(hspace=0.3,wspace=0.3,top=0.95,bottom=0.08,left=0.08,right=0.95)
    if not plot:
        plt.close()

    table2 = table2.sort_values(by='HZ_mp_min_osc+gr_texp15')
    print(table2[0:30][['ra_j2000','dec_j2000','primary_name','gmag','eff_nights','dist','HZ_mp_min_osc+gr_texp15']])
    
    return table2

"""
ra_j2000: --- Right ascension coordinate
dec_j2000: --- Declination coordinate
gmag: --- Gaia G magnitude
dist: --- Distance in parsec
Teff(K): --- Effective temperature of the star
logg: --- Surface gravity in cgs
Fe/H: --- Stellar mettalicity 
vsini: --- projected equatorial velocity in kms
airmass_min: --- minimum airmass along the year
tyr_rise: --- Decimal year start of the season
tyr_set: --- Decimal year end of the season
season_length: --- geometrical season length at La Palma above airmass of 1.75
eff_nights: --- Effective number of nights (Product of weather forecast with geometrical season length)
eff_airmass_mean: --- Mean optimal airmass value (Product of weather forecast with visibility curve)
eff_seeing: --- Mean expected seeing value (Function of mean_airmass and seasonal seeing)
nb_subexp: --- Number of subexp required to not saturate the detector on 15 minutes exposures
snr_550: 
snr_550_texp550: --- Expected SNR at 550 nm for 15min exposure (Clark ETC + Atmos parameters)
sig_rv_texp15: --- Expected RV precision at 550 nm for 15min exposure (Clark ETC + Atmos parameters)
texp_snr_250: --- Required exposure time to reach a SNR=250 at 550 nm
sig_rv_osc_texp15: --- Expected RV residual due to p-mode with 15min exposure
sig_rv_osc+gr_texp15: --- Expected RV residual due to p-mode + granulation with 15min exposure
HZ_period_inf: --- Inner habitable zone period
HZ_period_inf: --- outer habitable zone period
HZ_amp_inf: --- Inner habitable zone RV K-semi amplitude (1 earth mass)
HZ_amp_sup: --- Outer habitable zone RV K-semi amplitude (1 earth mass)
HZ_mp_min_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night for 8 years (photo noise)
HZ_mp_min_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night for 8 years (photo noise)
HZ_mp_min_osc_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night for 8 years (photo noise + p-mode)
HZ_mp_min_osc+gr_texp15: --- the minimum mass planet that can be detected in the inner HZ if observed every night (photo noise + p-mode + granulation)
"""


cutoff0 = {
    'Teff(K)<':6000,
    'Teff(K)>':4000,
    'logg>':4.0,
    'vsini<':8,
    'Fe/H>':-0.4,
    'Age>':9,
    'airmass_min>':0.99,
    'eff_nights>':160,
    'season_length>':180,
    'gmag<':7,
    'sig_rv_texp15>':0,
    'sig_rv_osc_texp15>':0,
    'sig_rv_osc+gr_texp15>':0,
    'HZ_mp_min_texp15>':0,
    'HZ_mp_min_osc+gr_texp15<':10,
    }

table_filtered = func_cutoff(table,cutoff0)

