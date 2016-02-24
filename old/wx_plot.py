# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

station_id = 'KCABAKER38'
full = 'COCONUT BAROMETER'
obs_comb = pd.read_pickle('/Users/areed145/Dropbox/wx_scraper/obs_data.pk1')

obs_comb['obs_d'] = pd.DatetimeIndex(obs_comb['obs_time']).tz_localize('UTC').tz_convert('US/Pacific').format(formatter=(lambda x: x.strftime('%Y%m%d')))
obs_comb['obs_time'] = pd.DatetimeIndex(obs_comb['obs_time']).tz_localize('UTC').tz_convert('US/Pacific')

obs_comb = obs_comb.convert_objects(convert_numeric=True)
obs_comb.loc[obs_comb.pressure_in < 10, 'pressure_in'] = pd.np.nan

obs_plot = obs_comb.set_index(['obs_time'])

obs_plot = obs_plot[-4500:]

obs_plot = obs_plot.resample('2min', how = 'mean').reset_index()
obs_plot = obs_plot.interpolate()

width = 2.5/(len(obs_plot))
    
fig = mpl.pyplot.gcf()
fig.set_size_inches(10.5, 16.5)

plt.subplot(8, 1, 1)
plt.title(station_id+': '+full)
plt.plot(obs_plot['obs_time'], obs_plot['temp_f'], 'b-')
plt.plot(obs_plot['obs_time'], obs_plot['dewpoint_f'], 'g-')

try:
    plt.plot(obs_plot['obs_time'], obs_plot['heat_index_f'], 'r-')
    plt.plot(obs_plot['obs_time'], obs_plot['wind_chill_f'], 'c-')
except:
    pass

plt.ylabel('Temp/Dewpt (F)')
plt.grid(b=True, which='both', color='k',linestyle='-')

plt.subplot(8, 1, 2)
plt.plot(obs_plot['obs_time'], obs_plot['pressure_in'], 'k-')
plt.ylabel('Pres (in)')
#plt.ylim([29.5,30.5])
plt.grid(b=True, which='both', color='k',linestyle='-')

plt.subplot(8, 1, 3)
plt.plot(obs_plot['obs_time'], obs_plot['wind_mph'], 'y-')
plt.plot(obs_plot['obs_time'], obs_plot['wind_gust_mph'], 'r.')
plt.ylabel('Wind/Gust (mph)')
plt.grid(b=True, which='both', color='k',linestyle='-')

obs_windNA = obs_plot[obs_plot.wind_mph == 0]
obs_wind00 = obs_plot[(obs_plot.wind_mph > 0) & (obs_plot.wind_mph <= 2)]
obs_wind02 = obs_plot[(obs_plot.wind_mph > 2) & (obs_plot.wind_mph <= 5)]
obs_wind05 = obs_plot[obs_plot.wind_mph >= 5]

obs_gustNA = obs_plot[obs_plot.wind_gust_mph == 0]
obs_gust00 = obs_plot[(obs_plot.wind_gust_mph > 0) & (obs_plot.wind_gust_mph <= 2)]
obs_gust02 = obs_plot[(obs_plot.wind_gust_mph > 2) & (obs_plot.wind_gust_mph <= 5)]
obs_gust05 = obs_plot[obs_plot.wind_gust_mph >= 5]

plt.subplot(8, 1, 4)
plt.bar(obs_gustNA['obs_time'], obs_gustNA['wind_degrees']-obs_gustNA['wind_degrees']+359, width, edgecolor='none', color='#FFFFFF')
plt.bar(obs_gust00['obs_time'], obs_gust00['wind_degrees']-obs_gust00['wind_degrees']+359, width, edgecolor='none', color='#FFCC00')
plt.bar(obs_gust02['obs_time'], obs_gust02['wind_degrees']-obs_gust02['wind_degrees']+359, width, edgecolor='none', color='#FF9900')
plt.bar(obs_gust05['obs_time'], obs_gust05['wind_degrees']-obs_gust05['wind_degrees']+359, width, edgecolor='none', color='#FF0000')
plt.bar(obs_windNA['obs_time'], obs_windNA['wind_degrees']-obs_windNA['wind_degrees']+250, width, edgecolor='none', color='#FFFFFF')
plt.bar(obs_wind00['obs_time'], obs_wind00['wind_degrees']-obs_wind00['wind_degrees']+250, width, edgecolor='none', color='#FFCC00')
plt.bar(obs_wind02['obs_time'], obs_wind02['wind_degrees']-obs_wind02['wind_degrees']+250, width, edgecolor='none', color='#FF9900')
plt.bar(obs_wind05['obs_time'], obs_wind05['wind_degrees']-obs_wind05['wind_degrees']+250, width, edgecolor='none', color='#FF0000')
plt.plot(obs_plot['obs_time'], obs_plot['wind_degrees'], 'k-', ms=2)
#plt.plot(obs_windNA['obs_time'], obs_windNA['wind_degrees'], '.', ms=3, color='#9191B5')
#plt.plot(obs_wind00['obs_time'], obs_wind00['wind_degrees'], '.', ms=5, color='#FFCC00')
#plt.plot(obs_wind02['obs_time'], obs_wind02['wind_degrees'], '.', ms=6, color='#FF9900')
#plt.plot(obs_wind05['obs_time'], obs_wind05['wind_degrees'], '.', ms=6, color='#FF0000')
plt.ylabel('Wind Dir (deg)')
plt.grid(b=True, which='both', color='k',linestyle='-')
plt.ylim([0,359])

plt.subplot(8, 1, 5)
plt.bar(obs_plot['obs_time'], obs_plot['cloud_base'], width, color='#99CCFF', edgecolor='none')
plt.ylabel('Min Cloudbase (ft MSL)')
plt.grid(b=True, which='both', color='k',linestyle='-')

plt.subplot(8, 1, 6)
plt.bar(obs_plot['obs_time'], obs_plot['solar_radiation']/100, width, color='#FFDB4D', edgecolor='none')
plt.bar(obs_plot['obs_time'], obs_plot['uv'], width, color='#DB4DB8', edgecolor='none')
plt.ylabel('Solar Rad (w/mm^2) / UV')
plt.grid(b=True, which='both', color='k',linestyle='-')

plt.subplot(8, 1, 7)
plt.bar(obs_plot['obs_time'], obs_plot['relative_humidity'], width, color='#256EB8', edgecolor='none')
plt.bar(obs_plot['obs_time'], 100 - obs_plot['relative_humidity'], width, bottom=obs_plot['relative_humidity'], color='#B2CCFF', edgecolor='none')
plt.ylabel('Rel Humidity (%)')
plt.grid(b=True, which='both', color='k',linestyle='-')

plt.subplot(8, 1, 8)
plt.bar(obs_plot['obs_time'], obs_plot['precip_today_in'], width, color='#0033CC', edgecolor='none')
plt.plot(obs_plot['obs_time'], obs_plot['precip_1hr_in'], 'c-')
plt.ylabel('Precip (in)')
plt.grid(b=True, which='both', color='k',linestyle='-')

plt.show()

fig.savefig('/Users/areed145/Dropbox/wx_scraper/obs_data_test.png', dpi=200)

wind_dir_mph = obs_comb.groupby(['obs_time', 'wind_dir']).wind_mph.mean().reset_index()