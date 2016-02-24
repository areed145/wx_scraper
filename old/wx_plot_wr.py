# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

station_id = 'KCABAKER38'
full = 'COBA'
cur_obs_date = '20151008'

obs_comb = pd.read_pickle('/Users/areed145/Dropbox/wx_scraper/old/obs_data.pk1 copy')
obs_new = pd.read_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/pickle_20151018.pk1')

obs_windNA = obs_new[obs_new.wind_mph == 0]
obs_wind00 = obs_new[(obs_new.wind_mph > 0) & (obs_new.wind_mph <= 2)]
obs_wind02 = obs_new[(obs_new.wind_mph > 2) & (obs_new.wind_mph <= 5)]
obs_wind05 = obs_new[(obs_new.wind_mph > 5) & (obs_new.wind_mph <= 10)]
obs_wind10 = obs_new[obs_new.wind_mph >= 10]

windNA = obs_windNA.groupby(['wind_dir'])['obs_d'].count().reset_index().rename(columns={'obs_d': 'windNA'})
wind00 = obs_wind00.groupby(['wind_dir'])['obs_d'].count().reset_index().rename(columns={'obs_d': 'wind00'})
wind02 = obs_wind02.groupby(['wind_dir'])['obs_d'].count().reset_index().rename(columns={'obs_d': 'wind02'})
wind05 = obs_wind05.groupby(['wind_dir'])['obs_d'].count().reset_index().rename(columns={'obs_d': 'wind05'})
wind10 = obs_wind10.groupby(['wind_dir'])['obs_d'].count().reset_index().rename(columns={'obs_d': 'wind10'})

df = windNA.merge(wind00, on=['wind_dir'], how='outer')\
            .merge(wind02, on=['wind_dir'], how='outer')\
            .merge(wind05, on=['wind_dir'], how='outer')\
            .merge(wind10, on=['wind_dir'], how='outer')
            
df = df.fillna(0)
order = {'North':0, 'NNE':1, 'NE':2, 'ENE':3, 'East':4, 'ESE':5, 'SE':6, 'SSE':7, 'South':8, 'SSW':9, 'SW':10, 'WSW':11, 'West':12, 'WNW':13, 'NW':14, 'NNW':15}  
df['order'] = df['wind_dir'].map(order)
df = df.sort(columns=['order'], axis=0, ascending=True)
            
total = df.windNA.sum()+df.wind00.sum()+df.wind02.sum()+df.wind05.sum()+df.wind10.sum()

df['calm'] = df.windNA.sum()/total*100/16
df.wind00 = df.wind00/total*100
df.wind02 = df.wind02/total*100
df.wind05 = df.wind05/total*100
df.wind10 = df.wind10/total*100
            
N = 16
theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
radii = 10 * np.random.rand(N)
width = 2 * np.pi / 16
theta = theta - (width/2)

fig = mpl.pyplot.gcf()
fig.set_size_inches(8, 8)
plt.suptitle(station_id+': '+full+'\n'+cur_obs_date, fontsize=10)
        
ax = plt.subplot(111, polar=True)
bars = ax.bar(theta, df.calm, width=width, color='#8AB8E6', edgecolor='none', bottom=0.0)
bars = ax.bar(theta, df.wind00, width=width, color='#FFCC00', edgecolor='#FFFFFF', bottom=df.calm)
bars = ax.bar(theta, df.wind02, width=width, color='#FF9900', edgecolor='#FFFFFF', bottom=df.calm+df.wind00)
bars = ax.bar(theta, df.wind05, width=width, color='#FF0000', edgecolor='#FFFFFF', bottom=df.calm+df.wind00+df.wind02)
bars = ax.bar(theta, df.wind10, width=width, color='#FF3385', edgecolor='#FFFFFF', bottom=df.calm+df.wind00+df.wind02+df.wind05)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)

## Use custom colors and opacity
#for r, bar in zip(radii, bars):
#    bar.set_alpha(0.5)

plt.show()
