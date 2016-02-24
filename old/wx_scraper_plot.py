# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import bs4
import urllib
import numpy as np
import pandas as pd
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
from cmath import rect, phase
from math import radians, degrees
from matplotlib.dates import DateFormatter
import pytz
import glob

sid = 'KCABAKER38'
resample_int = '5min'
width = 1/24/60*5 #interval is a day, need to divide by reample int
p_int = 32
folder = '/Users/areed145/Dropbox/wx_scraper/'
cur_obs_date = '20151019'
url = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID='+sid

try:
    soup = bs4.BeautifulSoup(urllib.request.urlopen(url))
    station_id = soup.find('station_id').getText()
    full = soup.find('full').getText()
except:
    pass

pickle = folder+'pickles/pickle_'+cur_obs_date+'.pk1'
archive = folder+'archives/archive_'+cur_obs_date+'.csv'
plot_obs = folder+'plots/obs_'+cur_obs_date+'.png'
plot_obs_cur = folder+'plots/obs_cur.png'

df_store = pd.DataFrame()

for pk1 in glob.glob(folder+'pickles/pickle_*.pk1'):
    df_store = df_store.append(pd.read_pickle(pk1))

last_obs_time = df_store.index.max()

def mean_angle(deg):
    if len(deg) > 2:
        deg = deg.astype(float)
        deg_avg = degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))
        if deg_avg < 0:
            deg_avg = 360 + deg_avg
        return deg_avg
    return pd.np.nan


obs_plot = df_store.resample(resample_int, how = 'mean')
obs_plot_wind = df_store[['wind_degrees']].resample(resample_int, how = mean_angle)
obs_plot['wind_degrees'] = obs_plot_wind['wind_degrees']
obs_plot_max = df_store[['wind_gust_mph']].resample(resample_int, how = 'max')
obs_plot['wind_gust_mph'] = obs_plot_max['wind_gust_mph']

obs_plot['wind_degrees_u'] = np.round((obs_plot['wind_degrees'] / (360 / p_int)),0) * (360 / p_int)
obs_plot.loc[obs_plot.wind_degrees_u == 360,'wind_degrees_u'] = 0
obs_windNA = obs_plot[obs_plot.wind_mph == 0]
obs_wind00 = obs_plot[(obs_plot.wind_mph > 0) & (obs_plot.wind_mph <= 1)]
obs_wind01 = obs_plot[(obs_plot.wind_mph > 1) & (obs_plot.wind_mph <= 2)]
obs_wind02 = obs_plot[(obs_plot.wind_mph > 2) & (obs_plot.wind_mph <= 5)]
obs_wind05 = obs_plot[(obs_plot.wind_mph > 5) & (obs_plot.wind_mph <= 10)]
obs_wind10 = obs_plot[obs_plot.wind_mph >= 10]
obs_gustNA = obs_plot[obs_plot.wind_gust_mph == 0]
obs_gust00 = obs_plot[(obs_plot.wind_gust_mph > 0) & (obs_plot.wind_gust_mph <= 1)]
obs_gust01 = obs_plot[(obs_plot.wind_gust_mph > 1) & (obs_plot.wind_gust_mph <= 2)]
obs_gust02 = obs_plot[(obs_plot.wind_gust_mph > 2) & (obs_plot.wind_gust_mph <= 5)]
obs_gust05 = obs_plot[(obs_plot.wind_gust_mph > 5) & (obs_plot.wind_gust_mph <= 10)]
obs_gust10 = obs_plot[obs_plot.wind_gust_mph >= 10]

plt.close("all")

fig = mpl.pyplot.gcf()
fig.set_size_inches(10, 30)
formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))

plt.rc('grid', color='1', linewidth=2, linestyle='-')

plt.suptitle(station_id+': '+full+'\n'+cur_obs_date+', Last Observation: '+last_obs_time.strftime('%H:%M:%S'), fontsize=20)

ax0 = plt.subplot2grid((9,2), (5, 0), colspan=2)
ax0.set_axis_bgcolor('.9')
ax0.set_axisbelow(True)
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.xaxis.set_ticks_position('none')
ax0.yaxis.set_ticks_position('none')
plt.plot_date(obs_plot.index, obs_plot['cloud_base']/100, marker = '', color='#99CCFF',linestyle='-', linewidth=3.0)
plt.ylabel('Min Cloudbase (ft MSL)')
#plt.ylim([0,100])
plt.grid(b=True, which='both', color='1',linestyle='-')

ax1 = plt.subplot2grid((9,2), (0, 0), colspan=2)
ax1.set_axis_bgcolor('.9')
ax1.set_axisbelow(True)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.xaxis.set_ticks_position('none')
ax1.yaxis.set_ticks_position('none')
plt.plot_date(obs_plot.index, obs_plot['temp_f'], marker = '', color='#3366FF',linestyle='-', linewidth=3.0)
plt.plot_date(obs_plot.index, obs_plot['dewpoint_f'], marker = '', color='#00CC00',linestyle='-', linewidth=3.0)

try:
    plt.plot_date(obs_plot.index, obs_plot['heat_index_f'], marker = '', color='#FF0066',linestyle='-', linewidth=3.0)
    plt.plot_date(obs_plot.index, obs_plot['wind_chill_f'], marker = '', color='#00FFCC',linestyle='-', linewidth=3.0)
except:
    pass

plt.ylabel('Temp/Dewpt (F)')
#plt.ylim([0,110])
plt.grid(b=True, which='both', color='1',linestyle='-')

ax2 = plt.subplot2grid((9,2), (1, 0), colspan=2)
ax2.set_axis_bgcolor('.9')
ax2.set_axisbelow(True)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.xaxis.set_ticks_position('none')
ax2.yaxis.set_ticks_position('none')
ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
plt.plot_date(obs_plot.index, obs_plot['pressure_in'], marker = '', color='#FF5C33',linestyle='-', linewidth=3.0)
plt.ylabel('Pres (in)')
#plt.ylim([29.6,30.6])
plt.grid(b=True, which='both', color='1',linestyle='-')

ax3 = plt.subplot2grid((9,2), (2, 0), colspan=2)
ax3.set_axis_bgcolor('.9')
ax3.set_axisbelow(True)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.xaxis.set_ticks_position('none')
ax3.yaxis.set_ticks_position('none')
plt.bar(obs_gustNA.index, obs_gustNA['wind_degrees']-obs_gustNA['wind_degrees']-80, width, edgecolor='none', color='#8AB8E6')
plt.bar(obs_gust00.index, obs_gust00['wind_degrees']-obs_gust00['wind_degrees']-80, width, edgecolor='none', color='#FFE066')
plt.bar(obs_gust01.index, obs_gust01['wind_degrees']-obs_gust01['wind_degrees']-80, width, edgecolor='none', color='#FFCC00')
plt.bar(obs_gust02.index, obs_gust02['wind_degrees']-obs_gust02['wind_degrees']-80, width, edgecolor='none', color='#FF9900')
plt.bar(obs_gust05.index, obs_gust05['wind_degrees']-obs_gust05['wind_degrees']-80, width, edgecolor='none', color='#FF0000')
plt.bar(obs_gust10.index, obs_gust10['wind_degrees']-obs_gust10['wind_degrees']-80, width, edgecolor='none', color='#FF3385')
plt.bar(obs_windNA.index, obs_windNA['wind_degrees']-obs_windNA['wind_degrees']-40, width, edgecolor='none', color='#8AB8E6')
plt.bar(obs_wind00.index, obs_wind00['wind_degrees']-obs_wind00['wind_degrees']-40, width, edgecolor='none', color='#FFE066')
plt.bar(obs_wind01.index, obs_wind01['wind_degrees']-obs_wind01['wind_degrees']-40, width, edgecolor='none', color='#FFCC00')
plt.bar(obs_wind02.index, obs_wind02['wind_degrees']-obs_wind02['wind_degrees']-40, width, edgecolor='none', color='#FF9900')
plt.bar(obs_wind05.index, obs_wind05['wind_degrees']-obs_wind05['wind_degrees']-40, width, edgecolor='none', color='#FF0000')
plt.bar(obs_wind10.index, obs_wind10['wind_degrees']-obs_wind10['wind_degrees']-40, width, edgecolor='none', color='#FF3385')
plt.plot_date(obs_plot.index, obs_plot['wind_degrees'], 'b.', ms=7)
plt.yticks(np.linspace(0,360,(p_int/4)+1))
plt.ylabel('Wind Dir (deg)')
plt.ylim([-80,360])
plt.grid(b=True, which='both', color='1',linestyle='-')

ax4 = plt.subplot2grid((9,2), (3, 0), colspan=2)
ax4.set_axis_bgcolor('.9')
ax4.set_axisbelow(True)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.xaxis.set_ticks_position('none')
ax4.yaxis.set_ticks_position('none')
plt.bar(obs_plot.index, obs_plot['solar_radiation']/100, width, color='#FFDB4D', edgecolor='none')
plt.plot_date(obs_plot.index, obs_plot['uv'], marker = '', color='#DB4DB8',linestyle='-', linewidth=3.0)
plt.ylabel('Solar Rad (w/mm^2) / UV')
plt.ylim([0,10])
plt.grid(b=True, which='both', color='1',linestyle='-')

ax5 = plt.subplot2grid((9,2), (4, 0), colspan=2)
ax5.set_axis_bgcolor('.9')
ax5.set_axisbelow(True)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.spines['bottom'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax5.xaxis.set_ticks_position('none')
ax5.yaxis.set_ticks_position('none')
plt.plot_date(obs_plot.index, obs_plot['relative_humidity'], marker = '', color='#256EB8',linestyle='-', linewidth=3.0)
plt.ylabel('Rel Humidity (%)')
plt.ylim([0,100])
plt.grid(b=True, which='both', color='1',linestyle='-')

ax6 = plt.subplot2grid((9,2), (6, 0), colspan=2)
ax6.set_axis_bgcolor('.9')
ax6.set_axisbelow(True)
ax6.spines['top'].set_visible(False)
ax6.spines['right'].set_visible(False)
ax6.spines['bottom'].set_visible(False)
ax6.spines['left'].set_visible(False)
ax6.xaxis.set_ticks_position('none')
ax6.yaxis.set_ticks_position('none')
plt.bar(obs_plot.index, obs_plot['precip_today_in'], width, color='#335CD6', edgecolor='none')
plt.plot_date(obs_plot.index, obs_plot['precip_1hr_in'], marker = '', color='#00CCFF',linestyle='-', linewidth=3.0)
plt.ylabel('Precip (in)')
plt.ylim([0,0.5])
plt.grid(b=True, which='both', color='1',linestyle='-')

windNA = obs_windNA.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'windNA'})
wind00 = obs_wind00.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind00'})
wind01 = obs_wind01.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind01'})
wind02 = obs_wind02.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind02'})
wind05 = obs_wind05.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind05'})
wind10 = obs_wind10.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind10'})

df_deg = pd.DataFrame(data = np.linspace(0, 360, p_int, endpoint=False), columns = ['wind_degrees_u'])

df = df_deg.merge(windNA, on=['wind_degrees_u'], how='left')\
            .merge(wind00, on=['wind_degrees_u'], how='left')\
            .merge(wind01, on=['wind_degrees_u'], how='left')\
            .merge(wind02, on=['wind_degrees_u'], how='left')\
            .merge(wind05, on=['wind_degrees_u'], how='left')\
            .merge(wind10, on=['wind_degrees_u'], how='left')

df = df.fillna(0)
df = df.sort(columns=['wind_degrees_u'], axis=0, ascending=True)

total = df.windNA.sum()+df.wind00.sum()+df.wind01.sum()+df.wind02.sum()+df.wind05.sum()+df.wind10.sum()

df['calm'] = df.windNA.sum()/total*100/16
df.wind00 = df.wind00/total*100
df.wind01 = df.wind01/total*100
df.wind02 = df.wind02/total*100
df.wind05 = df.wind05/total*100
df.wind10 = df.wind10/total*100

theta = np.linspace(0.0, 2 * np.pi, p_int, endpoint=False)
width_polar = 2 * np.pi / p_int
theta = theta - (width_polar/2)

ax7 = plt.subplot2grid((9,2), (7, 0), colspan=2, rowspan=2, projection="polar")
ax7.set_axis_bgcolor('.9')
ax7.set_axisbelow(True)
ax7.spines['polar'].set_visible(False)
ax7.xaxis.set_ticks_position('none')
ax7.yaxis.set_ticks_position('none')
bars = ax7.bar(theta, df.calm, width=width_polar, color='#8AB8E6', edgecolor='none', bottom=0.0)
bars = ax7.bar(theta, df.wind00, width=width_polar, color='#FFE066', edgecolor='.9', bottom=df.calm)
bars = ax7.bar(theta, df.wind01, width=width_polar, color='#FFCC00', edgecolor='.9', bottom=df.calm+df.wind00)
bars = ax7.bar(theta, df.wind02, width=width_polar, color='#FF9900', edgecolor='.9', bottom=df.calm+df.wind00+df.wind01)
bars = ax7.bar(theta, df.wind05, width=width_polar, color='#FF0000', edgecolor='.9', bottom=df.calm+df.wind00+df.wind01+df.wind02)
bars = ax7.bar(theta, df.wind10, width=width_polar, color='#FF3385', edgecolor='.9', bottom=df.calm+df.wind00+df.wind01+df.wind02+df.wind05)
ax7.set_theta_zero_location('N')
ax7.set_theta_direction(-1)

ax0.xaxis.set_major_formatter(formatter)
ax1.xaxis.set_major_formatter(formatter)
ax2.xaxis.set_major_formatter(formatter)
ax3.xaxis.set_major_formatter(formatter)
ax4.xaxis.set_major_formatter(formatter)
ax5.xaxis.set_major_formatter(formatter)
ax6.xaxis.set_major_formatter(formatter)

fig.savefig(plot_obs, dpi=200)