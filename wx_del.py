# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""

import pandas as pd

df = pd.read_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/obs_data.pk1')
df2 = pd.read_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/pickle_20151018.pk1')
df.obs_time = pd.DatetimeIndex(df.obs_time, tz="GMT").tz_convert('US/Pacific')
df = df.convert_objects(convert_numeric=True)
df.ix[df.pressure_in < 28, 'pressure_in'] = pd.np.nan
df['cloud_base'] = (((df['temp_f'] - df['dewpoint_f']) / 4.4) * 1000) + 400
df.obs_d = df.obs_time.apply(lambda x: x.strftime('%Y%m%d'))
df = df.set_index('obs_time')
            
df[df.obs_d == '20151016'].to_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/pickle_20151016.pk1')
df[df.obs_d == '20151017'].to_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/pickle_20151017.pk1')

#df = df.drop('wind_degrees_u', axis=1)
#df.to_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/obs_data.pk1')