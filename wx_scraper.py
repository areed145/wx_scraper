# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import urllib
import time
import os
import gc
import re
import fnmatch
from io import StringIO
import datetime as dt
from cmath import rect, phase
from math import radians, degrees
from sys import platform as _platform
import numpy as np
from numpy import pi
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator, DayLocator, HourLocator
from dateutil.relativedelta import relativedelta
import bs4
from dropbox.client import DropboxClient

# input data
START_DATE = '2016-02-26' # date to start pulling data
SID_LIST = ['KCABAKER38', 'KCASANTA278', 'KCODENVE86', 'KCAINYOK7',
            'KTXDALLA233', 'MAU562', 'KCABAKER8', 'KDCWASHI122',
            'KFLKEYWE22', 'KWASEATT398', 'KNMSANTA66'] # list of stations to pull
SID_LIST = ['KCABAKER38'] # list of stations to pull
MAC_FOLDER = '/Users/areed145/Dropbox/GitHub/wx_scraper/' # FOLDER if on Mac
LINUX_FOLDER = '/home/pi/Desktop/wx_scraper/' # FOLDER if on Linux
WIN_FOLDER = 'C:/Users/bvjs/Python/python-3.4.3.amd64/bvjs/wx_scraper/' # FOLDER if on PC
P_INT = 16 # number of segments for wind rose plots
HEIGHT = 90 # HEIGHT of wind plot heat map
SUMM_TDY = 15 # minutes to aggregate in "today" summary
SUMM_DAY = 15 # minutes to aggregate in "day" summary
SUMM_WEK = 30 # minutes to aggregate in "week" summary
SUMM_MON = 2 # hours to aggregate in "month" summary
SUMM_3MO = 6 # hours ro aggregate in "3month" summary
SUMM_ALL = 24 # hours to aggreagate in "all" summary
SAVE_ARCHIVE = False # save archive?
SLEEP_TIME = 2 # minutes to sleep before pulling stations again
DROPBOX = True

# plot variables
MARKER_SIZE = 75
LINE_WIDTH = 1

# dropbox
ACCESS_TOKEN = '9-dK9AG5fnkAAAAAAAAySl9sqq6pdHP84Fx5jhGn8adduEqGcmI_xvUS33fwhVO5'

# convert hour summary intervals to minutes
SUMM_MON = SUMM_MON * 60
SUMM_3MO = SUMM_3MO * 60
SUMM_ALL = SUMM_ALL * 60

# update plot defaults
plt.rcParams.update({'font.size': 9})
plt.rcParams.update({'savefig.dpi': 200})
plt.rcParams.update({'savefig.facecolor': 'azure'})
plt.rcParams.update({'savefig.edgecolor': 'k'})
plt.rcParams.update({'savefig.format': 'png'})
plt.rcParams.update({'savefig.jpeg_quality': 95})
plt.rcParams.update({'savefig.pad_inches': 0.05})
plt.rcParams.update({'text.color': 'k'})

# set FOLDER
if _platform == "darwin":
    FOLDER = MAC_FOLDER
elif _platform == "win32":
    FOLDER = WIN_FOLDER
elif _platform == 'linux':
    FOLDER = LINUX_FOLDER

#check directories
DATA_DIR = FOLDER+'data/'
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)
ARCHIVE_DIR = FOLDER+'archive/'
if not os.path.exists(ARCHIVE_DIR):
    os.makedirs(ARCHIVE_DIR)

# define functions
def dropbox_upload(timef):
    """uploads plots to dropbox"""
    client = DropboxClient(ACCESS_TOKEN)

    for root, dirs, files in os.walk(plots_dir):
        
        includes = ['*_'+timef+'.*']
        includes = r'|'.join([fnmatch.translate(x) for x in includes])
        files = [f for f in files if re.match(includes, f)]
    
        for filename in files:
            # construct the full local path
            local_path = os.path.join(root, filename)    
            # construct the full Dropbox path
            relative_path = os.path.relpath(local_path, plots_dir)
            dropbox_path = os.path.join('/plots/'+SID, relative_path)    
            # upload the file
            try:
                with open(local_path, 'rb') as f:
                    client.put_file(dropbox_path, f, overwrite=True)
            except:
                pass

def mean_angle(deg):
    """returns average angle for wind direction (since in radians)"""
    deg = deg[~deg.isnull()]
    if len(deg) > 2:
        deg = deg.astype(float)
        deg_avg = degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))
        if deg_avg < 0:
            deg_avg = 360 + deg_avg
        return deg_avg
    return np.mean(deg)
    
def summary_table(data, timef):
    """summarizes data into table form"""
    day_min = pd.DataFrame(data=data.min(), columns=['Min'])
    day_mean = pd.DataFrame(data=data.mean(), columns=['Avg'])
    day_max = pd.DataFrame(data=data.max(), columns=['Max'])
    data = pd.merge(day_min, pd.merge(day_mean, day_max, left_index=True, right_index=True),
                    left_index=True, right_index=True)
    return data

def summarize(data, begin, r_int):
    """summarizes data into specified time windows and stats by time interval"""
    data_mean = data[data.Time >= begin].resample(str(r_int)+'min', how='mean')
    data_mean.rename(columns=lambda x: x+'_avg', inplace=True)
    data_min = data[data.Time >= begin].resample(str(r_int)+'min', how='min')
    data_min.rename(columns=lambda x: x+'_min', inplace=True)
    data_max = data[data.Time >= begin].resample(str(r_int)+'min', how='max')
    data_max.rename(columns=lambda x: x+'_max', inplace=True)
    data_mean['WindDir'] = data['WindDir'].resample(str(r_int)+'min', how=mean_angle)
    data = pd.merge(data_mean, pd.merge(data_min, data_max, left_index=True, right_index=True),
                    left_index=True, right_index=True)
    data.drop(['DateUTC_max',
               'DateUTC_min',
               'Time_max',
               'Time_min',
               'WindDir_avg',
               'WindDir_max',
               'WindDir_min',
               'WindDirection_max',
               'PrecipDaily_min',
               'PrecipHourly_min',
               'WindDirection_min'], axis=1, inplace=True)
    data = data.reindex_axis(sorted(data.columns), axis=1)
    data['Precip_cum'] = data.PrecipHourly_avg.cumsum() * r_int / 60
    for col in ['Dewpoint_avg', 'Dewpoint_max', 'Dewpoint_min',
                'PrecipDaily_avg', 'PrecipDaily_max',
                'PrecipHourly_avg', 'PrecipHourly_max',
                'Humidity_avg','Humidity_max','Humidity_min',
                'CloudBase_avg','CloudBase_max','CloudBase_min',
                'SolarRadiation_avg', 'SolarRadiation_max', 'SolarRadiation_min',
                'Temp_avg', 'Temp_max', 'Temp_min',
                'WindDir',
                'WindspeedGust_avg', 'WindspeedGust_max', 'WindspeedGust_min',
                'Windspeed_avg', 'Windspeed_max', 'Windspeed_min',
                'Precip_cum']:    
        data[col] = np.round(data[col], 1)        
    for col in ['Pressure_avg', 'Pressure_max', 'Pressure_min']:
        data[col] = np.round(data[col], 2)
        
    return data

def rawlimit_date(data, begin):
    """limits raw data by the data provided"""
    data = data[data.Time >= begin]
    data = data.reindex_axis(sorted(data.columns), axis=1)
    return data

def rawlimit_daynite(data):
    """limits the raw data into night and day"""
    data_day = data[(data.Time.map(lambda x: x.hour) >= 6) &
                    (data.Time.map(lambda x: x.hour) < 18)]
    data_nit = data[~data.Time.isin(data_day.Time)]
    data_day = data_day.reindex_axis(sorted(data_day.columns), axis=1)
    data_nit = data_nit.reindex_axis(sorted(data_nit.columns), axis=1)
    return data_day, data_nit

def rawlimit_season(data):
    """limits the data into seasonl data"""
    data_spr = data[(data.Time.map(lambda x: x.month) >= 3) &
                    (data.Time.map(lambda x: x.month) <= 5)]
    data_smr = data[(data.Time.map(lambda x: x.month) >= 6) &
                    (data.Time.map(lambda x: x.month) <= 8)]
    data_fal = data[(data.Time.map(lambda x: x.month) >= 9) &
                    (data.Time.map(lambda x: x.month) <= 11)]
    data_not = data[(data.Time.map(lambda x: x.month) >= 3) &
                    (data.Time.map(lambda x: x.month) <= 11)]
    data_win = data[~data.Time.isin(data_not.Time)]
    data_win = data_win.reindex_axis(sorted(data_win.columns), axis=1)
    data_spr = data_spr.reindex_axis(sorted(data_spr.columns), axis=1)
    data_smr = data_smr.reindex_axis(sorted(data_smr.columns), axis=1)
    data_fal = data_fal.reindex_axis(sorted(data_fal.columns), axis=1)
    return data_win, data_spr, data_smr, data_fal

def treat(data, col1, valmin, val, col2):
    """treats data by replacing values"""
    try:
        data.loc[data[col1] <= valmin, col1] = val
        data.rename(columns={col1:col2}, inplace=True)
    except ValueError:
        pass
    
def export(data, ftype, name, timef):
    """collects and exports all the files meeting search criteria to dropbox"""
    if ftype == 'raw':
        data['Date'] = data.Time.map(lambda x: pd.to_datetime(x).date())
        data['TimeH'] = data.Time.map(lambda x: pd.to_datetime(x).time())
        data.columns = ['Cloud Base (ft)', 'Date UTC', 'Dewpoint (F)', 'Humidity (%)',
                        'Precip Daily (in)', 'Precip Hourly (in)', 'Pressure (inHg)',
                        'Solar Radiation (W/m^2)', 'Temp (F)', 'Date / Time', 'Wind Deg (deg)',
                        'Wind Dir', 'Wind Speed (mph)', 'Wind Speed Gust (mph)',
                        'Date', 'Time']
        data = data[['Date / Time','Date','Time', 'Date UTC', 'Temp (F)', 'Dewpoint (F)',
                     'Humidity (%)', 'Pressure (inHg)', 'Wind Dir', 'Wind Deg (deg)',
                     'Wind Speed (mph)', 'Wind Speed Gust (mph)', 'Solar Radiation (W/m^2)',
                     'Cloud Base (ft)', 'Precip Daily (in)', 'Precip Hourly (in)']]
        data = data.sort(['Date / Time'], ascending=False)
    elif ftype == 'summary':
        data = data.reset_index()
        data['Date'] = data.Time.map(lambda x: pd.to_datetime(x).date())
        data['TimeH'] = data.Time.map(lambda x: pd.to_datetime(x).time())
        data.columns = ['Date / Time',
                        'Avg Cloud Base','Max Cloud Base','Min Cloud Base',
                        'Avg Dewpoint','Max Dewpoint','Min Dewpoint',
                        'Avg Humidity','Max Humidity','Min Humidity',
                        'Avg Precip Daily','Max Precip Daily',
                        'Avg Precip Hourly','Max Precip Hourly',
                        'Avg Pressure','Max Pressure','Min Pressure',
                        'Avg Solar Radiation','Max Solar Radiation','Min Solar Radiation',
                        'Avg Temp','Max Temp','Min Temp',
                        'Wind Deg',
                        'Avg Wind Speed Gust','Max Wind Speed Gust','Min Wind Speed Gust',
                        'Avg Wind Speed','Max Wind Speed','Min Wind Speed',
                        'Cumulative Precip','Date','Time']
        data = data[['Date / Time','Date','Time',
                    'Avg Temp','Max Temp','Min Temp',
                    'Avg Dewpoint','Max Dewpoint','Min Dewpoint',
                    'Avg Humidity','Max Humidity','Min Humidity',
                    'Avg Pressure','Max Pressure','Min Pressure',
                    'Wind Deg',
                    'Avg Wind Speed Gust','Max Wind Speed Gust','Min Wind Speed Gust',
                    'Avg Wind Speed','Max Wind Speed','Min Wind Speed',
                    'Avg Solar Radiation','Max Solar Radiation','Min Solar Radiation',
                    'Avg Cloud Base','Max Cloud Base','Min Cloud Base',
                    'Avg Precip Daily','Max Precip Daily',
                    'Avg Precip Hourly','Max Precip Hourly',
                    'Cumulative Precip']]
        data = data.sort(['Date / Time'], ascending=False)


    data.to_csv(plots_dir+SID+'_'+name+'_'+timef+'.csv')            

def wind_rose(data, name):
    """creates wind rose diagrams"""
    data.loc[:, 'WindDir'] = np.round((data.WindDir / (360 / P_INT)), 0) * (360 / P_INT)
    data.loc[data['WindDir'] == 360, 'WindDir'] = 0
    windna = data[data.Windspeed == 0]
    wind00 = data[(data.Windspeed > 0) & (data.Windspeed <= 1)]
    wind01 = data[(data.Windspeed > 1) & (data.Windspeed <= 2)]
    wind02 = data[(data.Windspeed > 2) & (data.Windspeed <= 5)]
    wind05 = data[(data.Windspeed > 5) & (data.Windspeed <= 10)]
    wind10 = data[data.Windspeed >= 10]
    gb_wind00 = wind00.groupby(['WindDir'])['Windspeed']\
                .count().reset_index().rename(columns={'Windspeed': 'wind00'})
    gb_wind01 = wind01.groupby(['WindDir'])['Windspeed']\
                .count().reset_index().rename(columns={'Windspeed': 'wind01'})
    gb_wind02 = wind02.groupby(['WindDir'])['Windspeed']\
        .count().reset_index().rename(columns={'Windspeed': 'wind02'})
    gb_wind05 = wind05.groupby(['WindDir'])['Windspeed']\
                .count().reset_index().rename(columns={'Windspeed': 'wind05'})
    gb_wind10 = wind10.groupby(['WindDir'])['Windspeed']\
                .count().reset_index().rename(columns={'Windspeed': 'wind10'})
    data_deg = pd.DataFrame(data=np.linspace(0, 360, P_INT, endpoint=False),
                            columns=['WindDir'])
    data_wr = data_deg.merge(gb_wind00, on=['WindDir'], how='left')\
    .merge(gb_wind01, on=['WindDir'], how='left')\
    .merge(gb_wind02, on=['WindDir'], how='left')\
    .merge(gb_wind05, on=['WindDir'], how='left')\
    .merge(gb_wind10, on=['WindDir'], how='left')
    data_wr = data_wr.fillna(0)
    data_wr = data_wr.sort(columns=['WindDir'], axis=0, ascending=True)
    total = len(windna)+len(wind00)+len(wind01)+len(wind02)+len(wind05)+len(wind10)
    data_wr = data_wr.fillna(0)
    data_wr = data_wr.sort(columns=['WindDir'], axis=0, ascending=True)
    total = len(windna)+len(wind00)+len(wind01)+len(wind02)+len(wind05)+len(wind10)
    data_wr['windna'] = len(windna)/total*100/P_INT
    data_wr.wind00 = data_wr.wind00/total*100
    data_wr.wind01 = data_wr.wind01/total*100
    data_wr.wind02 = data_wr.wind02/total*100
    data_wr.wind05 = data_wr.wind05/total*100
    data_wr.wind10 = data_wr.wind10/total*100
    theta = np.linspace(0.0, 2 * pi, P_INT, endpoint=False)
    width_polar = 2 * pi / P_INT
    theta = theta - (width_polar/2)
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(4.5, 4.5)
    ax1 = plt.subplot(projection="polar")
    ax1.set_axisbelow(True)
    ax1.spines['polar'].set_visible(True)
    ax1.bar(theta,
            data_wr.windna,
            width=width_polar,
            color='#3366ff',
            edgecolor='none',
            bottom=0.0)
    ax1.bar(theta,
            data_wr.wind00,
            width=width_polar,
            color='#009999',
            edgecolor='k',
            bottom=data_wr.windna)
    ax1.bar(theta,
            data_wr.wind01,
            width=width_polar,
            color='#00cc00',
            edgecolor='k',
            bottom=data_wr.windna\
                   +data_wr.wind00)
    ax1.bar(theta,
            data_wr.wind02,
            width=width_polar,
            color='#bfff00',
            edgecolor='k',
            bottom=data_wr.windna\
                   +data_wr.wind00\
                   +data_wr.wind01)
    ax1.bar(theta,
            data_wr.wind05,
            width=width_polar,
            color='#ffcc00',
            edgecolor='k',
            bottom=data_wr.windna\
                   +data_wr.wind00\
                   +data_wr.wind01\
                   +data_wr.wind02)
    ax1.bar(theta,
            data_wr.wind10,
            width=width_polar,
            color='#ffff00',
            edgecolor='k',
            bottom=data_wr.windna\
                   +data_wr.wind00\
                   +data_wr.wind01\
                   +data_wr.wind02\
                   +data_wr.wind05)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)
    fig.text(0.02, 0.98,
             SID+' - '+city+', '+state+'\n'\
                 +str(lat)+', '+str(long)+', '+str(elev)+\
                 'ft\nWind Rose - '+name,
             verticalalignment='top',
             horizontalalignment='left',
             transform=ax1.transAxes)
    fig.text(0.98, 0.02,
             latest,
             fontsize=7,
             verticalalignment='bottom',
             horizontalalignment='right',
             transform=ax1.transAxes)
    fig.savefig(plots_dir+SID+'_wr_'+name+'.'+plt.rcParams['savefig.format'])
    plt.close(fig)

def wind_date(data, name, label, major, minor):
    """creates wind date plots"""
    width = 1/24/60*int(data.index.freqstr[:-1])
    data.loc[:, 'WindDir'] = np.round((data.WindDir / (360 / P_INT)), 0) * (360 / P_INT)
    data.loc[data['WindDir'] == 360, 'WindDir'] = 0
    windna = data[data.Windspeed_avg == 0]
    wind00 = data[(data.Windspeed_avg > 0) & (data.Windspeed_avg <= 1)]
    wind01 = data[(data.Windspeed_avg > 1) & (data.Windspeed_avg <= 2)]
    wind02 = data[(data.Windspeed_avg > 2) & (data.Windspeed_avg <= 5)]
    wind05 = data[(data.Windspeed_avg > 5) & (data.Windspeed_avg <= 10)]
    wind10 = data[data.Windspeed_avg >= 10]
    gustna = data[data.WindspeedGust_max == 0]
    gust00 = data[(data.WindspeedGust_max > 0) & (data.WindspeedGust_max <= 1)]
    gust01 = data[(data.WindspeedGust_max > 1) & (data.WindspeedGust_max <= 2)]
    gust02 = data[(data.WindspeedGust_max > 2) & (data.WindspeedGust_max <= 5)]
    gust05 = data[(data.WindspeedGust_max > 5) & (data.WindspeedGust_max <= 10)]
    gust10 = data[data.WindspeedGust_max >= 10]
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(SID+' - '+city+', '+state+': '+str(lat)+\
        ', '+str(long)+', '+str(elev)+\
        'ft\nWind Plot - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    ax1.bar(gustna.index,
            np.zeros(len(gustna))-HEIGHT,
            width,
            edgecolor='none',
            color='#3366ff')
    ax1.bar(gust00.index,
            np.zeros(len(gust00))-HEIGHT,
            width,
            edgecolor='none',
            color='#009999')
    ax1.bar(gust01.index,
            np.zeros(len(gust01))-HEIGHT,
            width,
            edgecolor='none',
            color='#00cc00')
    ax1.bar(gust02.index,
            np.zeros(len(gust02))-HEIGHT,
            width,
            edgecolor='none',
            color='#bfff00')
    ax1.bar(gust05.index,
            np.zeros(len(gust05))-HEIGHT,
            width,
            edgecolor='none',
            color='#ffcc00')
    ax1.bar(gust10.index,
            np.zeros(len(gust10))-HEIGHT,
            width,
            edgecolor='none',
            color='#ffff00')
    ax1.bar(windna.index,
            np.zeros(len(windna))-(HEIGHT/2),
            width,
            edgecolor='none',
            color='#3366ff')
    ax1.bar(wind00.index,
            np.zeros(len(wind00))-(HEIGHT/2),
            width,
            edgecolor='none',
            color='#009999')
    ax1.bar(wind01.index,
            np.zeros(len(wind01))-(HEIGHT/2),
            width,
            edgecolor='none',
            color='#00cc00')
    ax1.bar(wind02.index,
            np.zeros(len(wind02))-(HEIGHT/2),
            width,
            edgecolor='none',
            color='#bfff00')
    ax1.bar(wind05.index,
            np.zeros(len(wind05))-(HEIGHT/2),
            width,
            edgecolor='none',
            color='#ffcc00')
    ax1.bar(wind10.index,
            np.zeros(len(wind10))-(HEIGHT/2),
            width,
            edgecolor='none',
            color='#ffff00')
    ax2 = ax1.twinx()
    ax2.plot_date(data.index,
                  data.WindspeedGust_max,
                  marker='',
                  color='seagreen',
                  linestyle='-',
                  linewidth=1)
    ax2.plot_date(data.index,
                  data.Windspeed_avg,
                  marker='',
                  color='yellowgreen',
                  linestyle='-',
                  linewidth=1)
    ylim2max = np.round(((data.WindspeedGust_max.max()+5) / 5), 0) * 5
    ylim2min = -ylim2max/4
    ax2.set_yticks(np.linspace(0, ylim2max, (ylim2max/4)+1))
    ax2.set_ylabel('Windspeed / Gust (mph)')
    ax2.yaxis.label.set_color('yellowgreen')
    ax1.plot_date(data.index, data.WindDir, 'rx', ms=5*LINE_WIDTH)
    ax1.set_yticks(np.linspace(0, 360, (P_INT/4)+1))
    ax2.set_ylim([ylim2min, ylim2max])
    ax1.set_ylabel('Wind Direction (deg)')
    ax1.yaxis.label.set_color('red')
    ax1.set_ylim([-HEIGHT, 360])
    ax1.xaxis.set_major_locator(major)
    ax1.xaxis.set_major_formatter(label)
    ax1.xaxis.set_minor_locator(minor)
    ax1.grid(b=True, which='major', color='k', linestyle='-')
    fig.text(0.98, 0.02,
             latest,
             fontsize=7,
             verticalalignment='bottom',
             horizontalalignment='right',
             transform=ax1.transAxes)
    fig.savefig(plots_dir+SID+'_wd_'+name+'.'+plt.rcParams['savefig.format'])
    plt.close(fig)

def dtdt_solar_temp(data, name):
    """creates solar temp plots"""
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(SID+' - '+city+', '+state+': '+str(lat)+\
        ', '+str(long)+', '+str(elev)+\
        'ft\nSolar Radiation, dT/dt, Temp - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    axcb = ax1.scatter(data.SolarRadiation,
                       data.dTdt,
                       marker='x',
                       c=data.Temp,
                       vmin=-20,
                       vmax=120,
                       s=MARKER_SIZE,
                       cmap='jet',
                       alpha=0.8,
                       lw=1)
    ax1.set_xlabel('Solar Radiation (W/m^2)')
    ax1.set_ylabel('dT/dt (degF/hr)')
    ax1.grid(b=True, which='both', color='k', linestyle='-')
    cbar = plt.colorbar(axcb)
    cbar.set_label('Temp (degF)')
    fig.text(0.98, 0.02,
             latest,
             fontsize=7,
             verticalalignment='bottom',
             horizontalalignment='right',
             transform=ax1.transAxes)
    fig.savefig(plots_dir+SID+'_dTdts_'+name+'.'+plt.rcParams['savefig.format'])
    plt.close(fig)

def dtdt_date(data, name, label, major, minor):
    """creates dTdt date plots"""
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(SID+' - '+city+', '+state+': '+str(lat)+\
        ', '+str(long)+', '+str(elev)+\
        'ft\nTemp + dT/dt - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    vminmax = max(abs(data.dTdt_avg.max()), abs(data.dTdt_avg.min()))
    axcb = ax1.scatter(data.index,
                       data.Temp_avg,
                       c=data.dTdt_avg,
                       marker='x',
                       s=MARKER_SIZE,
                       cmap='jet',
                       vmin=-vminmax,
                       vmax=vminmax,
                       alpha=0.8,
                       lw=1)
    ax1.set_xlim([data.index.min(), data.index.max()])
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Temp (degF)')
    ax1.xaxis.set_major_locator(major)
    ax1.xaxis.set_major_formatter(label)
    ax1.xaxis.set_minor_locator(minor)
    ax1.grid(b=True, which='major', color='k', linestyle='-')
    cbar = plt.colorbar(axcb)
    cbar.set_label('dT/dt (degF/hr)')
    fig.text(0.98, 0.02,
             latest,
             fontsize=7,
             verticalalignment='bottom',
             horizontalalignment='right',
             transform=ax1.transAxes)
    fig.savefig(plots_dir+SID+'_dTdtd_'+name+'.'+plt.rcParams['savefig.format'])
    plt.close(fig)

def temp_dew_hum(data, name):
    """creates temp dewpt humidity plots"""
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(SID+' - '+city+', '+state+': '+str(lat)+\
        ', '+str(long)+', '+str(elev)+\
        'ft\nTemp, Dewpoint, Humidity - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    axcb = ax1.scatter(data.Dewpoint,
                       data.Temp,
                       c=data.Humidity,
                       marker='x',
                       s=MARKER_SIZE,
                       cmap='winter',
                       vmin=0,
                       vmax=100,
                       alpha=0.8,
                       lw=1)
    ax1.set_xlabel('Dewpoint (degF)')
    ax1.set_ylabel('Temp (degF)')
    ax1.grid(b=True, which='both', color='k', linestyle='-')
    cbar = plt.colorbar(axcb)
    cbar.set_label('Humidity (%)')
    fig.text(0.98, 0.02,
             latest,
             fontsize=7,
             verticalalignment='bottom',
             horizontalalignment='right',
             transform=ax1.transAxes)
    fig.savefig(plots_dir+SID+'_tdhs_'+name+'.'+plt.rcParams['savefig.format'])
    plt.close(fig)

def dpdt_wind_temp(data, name):
    """creates dPdt wind temp plots"""
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(SID+' - '+city+', '+state+': '+str(lat)+\
        ', '+str(long)+', '+str(elev)+\
        'ft\nSolar Radiation, dT/dt, Temp - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    axcb = ax1.scatter(data.Windspeed,
                       data.dPdt,
                       c=data.Temp,
                       marker='x',
                       vmin=-20,
                       vmax=120,
                       s=MARKER_SIZE,
                       cmap='jet',
                       alpha=0.8,
                       lw=1)
    ax1.set_xlabel('Windspeed (mph)')
    ax1.set_ylabel('dP/dt (inHg/hr)')
    ax1.grid(b=True, which='both', color='k', linestyle='-')
    cbar = plt.colorbar(axcb)
    cbar.set_label('Temp (degF)')
    fig.text(0.98, 0.02,
             latest,
             fontsize=7,
             verticalalignment='bottom',
             horizontalalignment='right',
             transform=ax1.transAxes)
    fig.savefig(plots_dir+SID+'_dPdts_'+name+'.'+plt.rcParams['savefig.format'])
    plt.close(fig)

def combo(data, name, label, major, minor):
    """creates combo plots"""
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 9)
    fig.suptitle(SID+' - '+city+', '+state+'\n'\
        +str(lat)+', '+str(long)+', '+str(elev)+\
        'ft\nCombo Plot - '+name)
    ax1 = plt.subplot(5, 1, 1)
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    ax1.plot_date(data.index,
                  data.Temp_avg,
                  marker='',
                  color='orangered',
                  linestyle='-',
                  linewidth=LINE_WIDTH)
    ax1.fill_between(data.index,
                     data.Temp_min,
                     data.Temp_max,
                     facecolor='orangered',
                     alpha=.2)
    ax1.set_ylabel('Temperature (degF)')
    ax1.yaxis.label.set_color('orangered')
    ax1.set_ylim([0, 120])
    ax12 = ax1.twinx()
    ax12.plot_date(data.index,
                   data.Dewpoint_avg,
                   marker='',
                   color='deepskyblue',
                   linestyle='-',
                   linewidth=LINE_WIDTH)
    ax12.fill_between(data.index,
                      data.Dewpoint_min,
                      data.Dewpoint_max,
                      facecolor='deepskyblue',
                      alpha=.2)
    ax12.set_ylabel('Dewpoint (degF)')
    ax12.yaxis.label.set_color('deepskyblue')
    ax12.set_ylim([0, 120])
    ax1.xaxis.set_major_locator(major)
    ax1.xaxis.set_major_formatter(label)
    ax1.xaxis.set_minor_locator(minor)
    ax1.grid(b=True, which='major', color='k', linestyle='-')
    ax2 = plt.subplot(5, 1, 2)
    ax2.plot_date(data.index,
                  data.PrecipHourly_max,
                  marker='',
                  color='c',
                  linestyle='-',
                  linewidth=LINE_WIDTH)
    ax2.fill_between(data.index,
                     data.PrecipHourly_max,
                     data.PrecipHourly_max - data.PrecipHourly_max,
                     facecolor='c',
                     alpha=.2)
    ax2.set_ylabel('Precip (in/hr)')
    ax2.yaxis.label.set_color('c')
    ax22 = ax2.twinx()
    ax22.plot_date(data.index,
                   data.Precip_cum,
                   marker='',
                   color='dodgerblue',
                   linestyle='-',
                   linewidth=LINE_WIDTH)
    ax22.set_ylabel('Precip (in)')
    ax22.yaxis.label.set_color('dodgerblue')
    max_precip_hr = data.PrecipHourly_max.max()
    if max_precip_hr > 0:
        ylim2max = max_precip_hr + 0.1
    else:
        ylim2max = 0.1
    ax2.set_ylim([0, ylim2max])
    max_precip_cm = data.Precip_cum.max()
    if max_precip_cm > 0:
        ylim22max = max_precip_cm + 0.1
    else:
        ylim22max = 0.1
    ax22.set_ylim([0, ylim22max])
    ax2.xaxis.set_major_locator(major)
    ax2.xaxis.set_major_formatter(label)
    ax2.xaxis.set_minor_locator(minor)
    ax2.grid(b=True, which='major', color='k', linestyle='-')
    ax3 = plt.subplot(5, 1, 3)
    ax3.plot_date(data.index,
                  data.Humidity_avg,
                  marker='',
                  color='gold',
                  linestyle='-',
                  linewidth=LINE_WIDTH)
    ax3.fill_between(data.index,
                     data.Humidity_min,
                     data.Humidity_max,
                     facecolor='gold',
                     alpha=.2)
    ax3.set_ylabel('Humidity (%)')
    ax3.yaxis.label.set_color('gold')
    ax3.set_ylim([0, 100])
    ax32 = ax3.twinx()
    ax32.plot_date(data.index,
                   data.Pressure_avg,
                   marker='',
                   color='forestgreen',
                   linestyle='-',
                   linewidth=LINE_WIDTH)
    ax32.set_ylabel('Pressure (inHg)')
    ax32.yaxis.label.set_color('forestgreen')
    ax3.xaxis.set_major_locator(major)
    ax3.xaxis.set_major_formatter(label)
    ax3.xaxis.set_minor_locator(minor)
    ax3.grid(b=True, which='major', color='k', linestyle='-')
    ax4 = plt.subplot(5, 1, 4)
    ax4.plot_date(data.index,
                  data.CloudBase_avg,
                  marker='',
                  color='deeppink',
                  linestyle='-',
                  linewidth=LINE_WIDTH)
    ax4.fill_between(data.index,
                     data.CloudBase_min,
                     data.CloudBase_max,
                     facecolor='deeppink',
                     alpha=.2)
    ax4.set_ylabel('Minimum Cloudbase (ft)')
    ax4.yaxis.label.set_color('deeppink')
    ax42 = ax4.twinx()
    try:
        ax42.plot_date(data.index,
                       data.SolarRadiation_avg,
                       marker='',
                       color='orange',
                       linestyle='-',
                       linewidth=LINE_WIDTH)
        ax42.fill_between(data.index,
                          data.SolarRadiation_min,
                          data.SolarRadiation_max,
                          facecolor='orange',
                          alpha=.2)
    except ValueError:
        pass
    ax42.set_ylabel('Solar Radiation (W/m^2)')
    ax42.yaxis.label.set_color('orange')
    ax4.xaxis.set_major_locator(major)
    ax4.xaxis.set_major_formatter(label)
    ax4.xaxis.set_minor_locator(minor)
    ax4.grid(b=True, which='major', color='k', linestyle='-')
    ax5 = plt.subplot(5, 1, 5)
    ax5.plot_date(data.index,
                  data.WindDir,
                  'rx',
                  ms=5*LINE_WIDTH)
    ax5.set_yticks(np.linspace(0, 360, (P_INT/4)+1))
    ax5.set_ylabel('Wind Direction (deg)')
    ax5.yaxis.label.set_color('red')
    ax52 = ax5.twinx()
    ax52.plot_date(data.index,
                   data.WindspeedGust_max,
                   marker='',
                   color='seagreen',
                   linestyle='-',
                   linewidth=LINE_WIDTH)
    ax52.plot_date(data.index,
                   data.Windspeed_avg,
                   marker='',
                   color='yellowgreen',
                   linestyle='-',
                   linewidth=LINE_WIDTH)
    ax52.set_ylabel('Windspeed / Gust (mph)')
    ax52.yaxis.label.set_color('yellowgreen')
    ax5.set_ylim([0, 360])
    ax5.xaxis.set_major_locator(major)
    ax5.xaxis.set_major_formatter(label)
    ax5.xaxis.set_minor_locator(minor)
    ax5.grid(b=True, which='major', color='k', linestyle='-')
    ax5.set_xlabel('Date')
    fig.text(0.98, 0.02,
             latest,
             fontsize=7,
             verticalalignment='bottom',
             horizontalalignment='right',
             transform=ax5.transAxes)
    fig.savefig(plots_dir+SID+'_combo_'+name+'.'+plt.rcParams['savefig.format'])
    plt.close(fig)

# main loop
SID = None
begin = dt.datetime.today()

while 1 < 2:
    clock_15min = begin+relativedelta(minutes=15)
    clock_1hr = begin+relativedelta(minutes=60)
    clock_6hrs = begin+relativedelta(hours=6)
                    
    try:
        for SID in SID_LIST:
            plots_dir = FOLDER+'plots/'+SID+'/'
            if not os.path.exists(plots_dir):
                os.makedirs(plots_dir)

            # get station data
            url_sd = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID='+SID
            soup_sd = bs4.BeautifulSoup(urllib.request.urlopen(url_sd))
            #full = soup_sd.find('full').getText()
            neigh = soup_sd.find('neighborhood').getText()
            city = soup_sd.find('city').getText()
            state = soup_sd.find('state').getText()
            lat = float(soup_sd.find('latitude').getText())
            long = float(soup_sd.find('longitude').getText())
            elev = int(soup_sd.find('elevation').getText()[:-3])
            latest = soup_sd.find('observation_time').getText()

            start = dt.datetime.strptime(START_DATE, "%Y-%m-%d").date()
            today = dt.datetime.today().date()
            lim_day = begin+relativedelta(days=-1)
            lim_fch = today+relativedelta(days=-2)
            lim_wek = today+relativedelta(days=-7)
            lim_mon = today+relativedelta(months=-1)
            lim_3mo = today+relativedelta(months=-3)
            lim_yer = today+relativedelta(years=-1)
            
            fetch_range = pd.DataFrame(columns=['Date'],
                                       data=pd.date_range(START_DATE, today).date)
            try:
                data = pd.read_pickle(DATA_DIR+SID+'_data.pk1')
                data_dates = pd.DataFrame(columns = ['Date'],
                                          data=data.Time.map(lambda x: pd.to_datetime(x).date()).unique())
                recent_range = pd.DataFrame(columns=['Date'],
                                            data=pd.date_range(lim_fch, today).date)
    
                fetch_dates = fetch_range[~fetch_range.Date.isin(data_dates.Date)]\
                                .append(recent_range).drop_duplicates()
            except Exception:
                data = pd.DataFrame()
                fetch_dates = fetch_range

            #print header data
            print('-- WX_SCRAPER --')
            print('(c) Adam Reeder')
            print('station id: '+SID)
            print('neighborhood: '+neigh)
            print('city: '+city)
            print('state: '+state)
            print('lat: '+str(lat))
            print('long: '+str(long))
            print('elevation: '+str(elev))
            print('start date: '+str(start))
            print('history length: '+str(len(fetch_range))+' days')
            print('fetch length: '+str(len(fetch_dates))+' days')
            print(latest)

            # data pull loop
            for date in fetch_dates.Date:
                year = date.year
                month = date.month
                day = date.day
                year_str = str(year)
                if month < 10:
                    month_str = '0'+str(month)
                else:
                    month_str = str(month)
                if day < 10:
                    day_str = '0'+str(day)
                else:
                    day_str = str(day)
                url = 'http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID='\
                    +SID+'&day='+str(day)+'&month='+str(month)+'&year='+str(year)+\
                    '&graphspan=day&format=1'
                soup = bs4.BeautifulSoup(urllib.request.urlopen(url)).text
                if soup.count('\n') > 1:
                    data = data.append(pd.read_csv(StringIO(soup), index_col=False))
                    print(year_str+'-'+month_str+'-'+day_str+': fetched')
                else:
                    print(year_str+'-'+month_str+'-'+day_str+': empty')
                    continue
            data = data.drop_duplicates(subset=['Time'])
            data.to_pickle(DATA_DIR+SID+'_data.pk1')
                    
            print('data files loaded')

            # compile to main dataframe
            data['Time'] = pd.to_datetime(data.Time)
            data['DateUTC'] = pd.to_datetime(data.DateUTC)
            data.index = data.Time
            data = data.sort(['DateUTC'], ascending=True)

            # data treatment
            treat(data, 'DewpointF', -50, pd.np.nan, 'Dewpoint')
            treat(data, 'HourlyPrecipIn', 0, 0, 'PrecipHourly')
            treat(data, 'Humidity', 0, pd.np.nan, 'Humidity')
            treat(data, 'PressureIn', 0, pd.np.nan, 'Pressure')
            treat(data, 'SolarRadiationWatts/m^2', -100, pd.np.nan, 'SolarRadiation')
            treat(data, 'TemperatureF', -100, pd.np.nan, 'Temp')
            treat(data, 'WindDirectionDegrees', -.01, pd.np.nan, 'WindDir')
            treat(data, 'WindSpeedGustMPH', -.01, pd.np.nan, 'WindspeedGust')
            treat(data, 'WindSpeedMPH', -.01, pd.np.nan, 'Windspeed')
            treat(data, 'dailyrainin', -.01, pd.np.nan, 'PrecipDaily')
            data.loc[data['Windspeed'] == 0, 'WindDir'] = pd.np.nan
            print('data treated')

            # calculations
            data['SID'] = SID
            data['Lat'] = lat
            data['Long'] = long
            data['Elev'] = elev
            data['City'] = city
            data['State'] = state
            data['CloudBase'] = np.round((((data['Temp'] - data['Dewpoint']) / 4.4) * 1000) + elev, 0)
            #data['dT'] = data.Temp - data.Temp.shift(1)
            #data['dP'] = data.Pressure - data.Pressure.shift(1)
            #data['dt'] = data.Time - data.Time.shift(1)
            #data['dt'] = data.dt.astype('timedelta64[s]')/(60*60)
            #data['dTdt'] = data['dT'] / data['dt']
            #data['dPdt'] = data['dP'] / data['dt']
            data.drop(['SoftwareType',
                       'City',
                       'Clouds',
                       'Conditions',
                       'SID',
                       'Lat',
                       'Long',
                       'Elev',
                       'State'], axis=1, inplace=True)
            print('calculations completed')

            # archive
            if SAVE_ARCHIVE:
                data.to_csv(ARCHIVE_DIR+SID+'_archive.csv')
                print('data archived')
            else:
                print('data archive not selected')

            # create x-axis date labels
            label_day = DateFormatter("%H:%M")
            label_wek = DateFormatter("%m/%d")
            label_mon = DateFormatter("%m/%d")
            label_3mo = DateFormatter("%m/%y")
            label_all = DateFormatter("%m/%y")
            maj_day = HourLocator(interval=3)
            maj_wek = DayLocator()
            maj_mon = DayLocator(interval=7)
            maj_3mo = MonthLocator()
            maj_all = MonthLocator()
            min_day = HourLocator()
            min_wek = HourLocator()
            min_mon = DayLocator()
            min_3mo = DayLocator()
            min_all = DayLocator()
            
            print('run began: '+str(begin))
            print('time now: '+str(dt.datetime.today()))
            print('time in 15min: '+str(clock_15min))
            print('time in 1hr: '+str(clock_1hr))
            print('time in 6hrs: '+str(clock_6hrs))      

            # create limited dataframes and summaries
            data_rawl_day = rawlimit_date(data, lim_day)
            data_rawl_tdy = rawlimit_date(data, today)
            wind_rose(data_rawl_day, '1day')
            #dtdt_solar_temp(data_rawl_day, '1day')
            temp_dew_hum(data_rawl_day, '1day')
            #dpdt_wind_temp(data_rawl_day, '1day')
            tdy_table = summary_table(data_rawl_tdy, 'Day').T
            day_table = summary_table(data_rawl_day, 'Day').T
            wek_table = summary_table(rawlimit_date(data, lim_wek), 'Week').T
            mon_table = summary_table(rawlimit_date(data, lim_mon), 'Month').T
            yer_table = summary_table(rawlimit_date(data, lim_yer), 'Year').T
            
            cols = ['Cloud Base (ft)', 'Dewpoint (F)', 'Humidity (%)',
                                'Precip Daily (in)', 'Precip Hourly (in)', 'Pressure (inHg)',
                                'Solar Radiation (W/m^2)', 'Temp (F)', 'Wind Deg (deg)',
                                'Wind Speed (mph)', 'Wind Speed Gust (mph)']
            
            tdy_table.columns = cols
            day_table.columns = cols
            wek_table.columns = cols
            mon_table.columns = cols
            yer_table.columns = cols
            
            tdy_table.apply(lambda x: x.astype(float).round(2)).to_csv(plots_dir+SID+'_tdy_table_1day.csv')
            day_table.apply(lambda x: x.astype(float).round(2)).to_csv(plots_dir+SID+'_day_table_1day.csv')
            wek_table.apply(lambda x: x.astype(float).round(2)).to_csv(plots_dir+SID+'_wek_table_1day.csv')
            mon_table.apply(lambda x: x.astype(float).round(2)).to_csv(plots_dir+SID+'_mon_table_1day.csv')
            yer_table.apply(lambda x: x.astype(float).round(2)).to_csv(plots_dir+SID+'_yer_table_1day.csv')

            print('summary table created')            
            
            export(data_rawl_day, 'raw', 'data_rawl_day', '1day')
            export(data_rawl_tdy, 'raw', 'data_rawl_tdy', '1day')
            del data_rawl_day, data_rawl_tdy
            
            data_summ_day = summarize(data, lim_day, SUMM_DAY)
            data_summ_tdy = summarize(data, today, SUMM_DAY)
            combo(data_summ_day, '1day', label_day, maj_day, min_day)
            wind_date(data_summ_day, '1day', label_day, maj_day, min_day)
            #dtdt_date(data_summ_day, '1day', label_day, maj_day, min_day)
            export(data_summ_day, 'summary', 'data_summ_day', '1day')
            export(data_summ_tdy, 'summary', 'data_summ_tdy', '1day')
            del data_summ_day, data_summ_tdy
            
            print('daily plots and summary table created')
            
            if DROPBOX == True:
                dropbox_upload('1day')
                    
                print('daily plots and summary table uploaded to dropbox')
            
            if dt.datetime.today() >= clock_15min:
                
                data_rawl_wek = rawlimit_date(data, lim_wek)
                wind_rose(data_rawl_wek, '1week')
                #dtdt_solar_temp(data_rawl_wek, '1week')
                #temp_dew_hum(data_rawl_wek, '1week')
                #dpdt_wind_temp(data_rawl_wek, '1week')
                #export(data_rawl_wek, 'raw', 'data_rawl_wek', '1week')
                del data_rawl_wek
                
                data_summ_wek = summarize(data, lim_wek, SUMM_WEK)
                combo(data_summ_wek, '1week', label_wek, maj_wek, min_wek)
                wind_date(data_summ_wek, '1week', label_wek, maj_wek, min_wek)
                #dtdt_date(data_summ_wek, '1week', label_wek, maj_wek, min_wek)
                export(data_summ_wek, 'summary', 'data_summ_wek', '1week')
                del data_summ_wek
                
                print('weekly plots created')
                
                if DROPBOX == True:
                    dropbox_upload('1week')
                    
                    print('weekly plots uploaded to dropbox')
                
            else:
                print('weekly plots skipped')
                
            if dt.datetime.today() >= clock_1hr:
    
                data_rawl_mon = rawlimit_date(data, lim_mon)
                wind_rose(data_rawl_mon, '1month')
                #dtdt_solar_temp(data_rawl_mon, '1month')
                #temp_dew_hum(data_rawl_mon, '1month')
                #dpdt_wind_temp(data_rawl_mon, '1month')
                #export(data_rawl_mon, 'raw', 'data_rawl_mon', '1month')
                del data_rawl_mon
                
                data_summ_mon = summarize(data, lim_mon, SUMM_MON)
                combo(data_summ_mon, '1month', label_mon, maj_mon, min_mon)
                wind_date(data_summ_mon, '1month', label_mon, maj_mon, min_mon)
                #dtdt_date(data_summ_mon, '1month')
                export(data_summ_mon, 'summary', 'data_summ_mon', '1month')
                del data_summ_mon
                
                print('monthly plots created')
                
                if DROPBOX == True:
                    dropbox_upload('1month')
                    
                    print('monthly plots uploaded to dropbox')
                
            else:
                print('monthly plots skipped')
                
            if dt.datetime.today() >= clock_6hrs:
    
                data_rawl_3mo = rawlimit_date(data, lim_3mo)
                wind_rose(data_rawl_3mo, '3month')
                #dtdt_solar_temp(data_rawl_3mo, '3month')
                #temp_dew_hum(data_rawl_3mo, '3month')
                #dpdt_wind_temp(data_rawl_3mo, '3month')
                #export(data_rawl_3mo, 'raw', 'data_rawl_3mo', '3month')
                del data_rawl_3mo
                
                data_summ_3mo = summarize(data, lim_3mo, SUMM_3MO)
                #combo(data_summ_3mo, '3month', label_3mo, maj_3mo, min_3mo)
                wind_date(data_summ_3mo, '3month', label_3mo, maj_3mo, min_3mo)
                #dtdt_date(data_summ_3mo, '3month')
                export(data_summ_3mo, 'summary', 'data_summ_3mo', '3month')
                del data_summ_3mo
                
                print('3-month plots created')
            
                if DROPBOX == True:
                    dropbox_upload('3month')
                    
                    print('3-month plots uploaded to dropbox')
                
            else:
                print('3-month plots skipped')
                
            if dt.datetime.today() >= clock_6hrs:
    
                data_rawl_yer = rawlimit_date(data, lim_yer)
                wind_rose(data_rawl_yer, '1year')
                #dtdt_solar_temp(data_rawl_yer, '1year')
                temp_dew_hum(data_rawl_yer, '1year')
                #dpdt_wind_temp(data_rawl_yer, '1year')
                #export(data_rawl_yer, 'raw', 'data_rawl_yer', '1year')
                del data_rawl_yer
                
                data_summ_yer = summarize(data, lim_yer, SUMM_ALL)
                combo(data_summ_yer, '1year', label_all, maj_all, min_all)
                wind_date(data_summ_yer, '1year', label_all, maj_all, min_all)
                #dtdt_date(data_summ_yer, '1year')
                export(data_summ_yer, 'summary', 'data_summ_yer', '1year')
                del data_summ_yer
                
                print('yearly plots created')
                
                if DROPBOX == True:
                    dropbox_upload('1year')
                    
                    print('yearly plots uploaded to dropbox')
                
            else:
                print('yearly plots skipped')
                
            if dt.datetime.today() >= clock_6hrs:
    
                data_rawl_all = rawlimit_date(data, START_DATE)
                wind_rose(data_rawl_all, 'all')
                #dtdt_solar_temp(data_rawl_all, 'all')
                #temp_dew_hum(data_rawl_all, 'all')
                #dpdt_wind_temp(data_rawl_all, 'all')
                del data_rawl_all
                
                #data_summ_all = summarize(data, START_DATE, SUMM_ALL)
                #combo(data_summ_all, 'all', label_all, maj_all, min_all)
                #wind_date(data_summ_all, 'all', label_all, maj_all, min_all)
                #dtdt_date(data_summ_all, 'all')
                #del data_summ_all
                
                data_rawl_dys, data_rawl_nts = rawlimit_daynite(data)
                wind_rose(data_rawl_dys, 'days')
                del data_rawl_dys
                wind_rose(data_rawl_nts, 'nights')
                del data_rawl_nts
                
                #data_rawl_win, data_rawl_spr, data_rawl_smr, data_rawl_fal = rawlimit_season(data)
                
                print('all-time plots created')
                
                if DROPBOX == True:
                    dropbox_upload('all')
                    dropbox_upload('days')
                    dropbox_upload('nights')
                    
                    print('all-time plots uploaded to dropbox')
                
                begin = dt.datetime.today()
                
                print('clock reset')
                
            else:
                print('all-time plots skipped')

            del data
            gc.collect
            
    except ValueError:
        print(SID+' failed')

    print('sleeping for '+str(SLEEP_TIME)+' minute(s)')
    time.sleep(60*SLEEP_TIME)
