# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import bs4
import urllib
import glob
import time
import os
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from cmath import rect, phase
from math import radians, degrees
from dateutil.relativedelta import relativedelta
from sys import platform as _platform
#from matplotlib.dates import DateFormatter
#import pytz

#####################################################################################################

# input data
start_date = '2015-02-26' # date to start pulling data
sid_list = ['KCABAKER38','KCABAKER8'] # list of stations to pull
mac_folder = '/Users/areed145/Dropbox/GitHub/wx_scraper/' # folder if on Mac
win_folder = 'C:/Users/bvjs/Python/python-3.4.3.amd64/bvjs/wx_scraper/' # folder if on PC
p_int = 16 # number of segments for wind rose plots
height = 90 # height of wind plot heat map
summ_tdy = 15 # minutes to aggregate in "today" summary
summ_wek = 30 # minutes to aggregate in "week" summary
summ_mon = 2 # hours to aggregate in "month" summary
summ_3mo = 6 # hours ro aggregate in "3month" summary
summ_all = 24 # hours to aggreagate in "all" summary
save_archive = False # save archive?
sleep_time = 10 # minutes to sleep before pulling stations again

# convert hour summary intervals to minutes
summ_mon = summ_mon * 60
summ_3mo = summ_3mo * 60
summ_all = summ_all * 60

# update plot defaults
plt.rcParams.update({'font.size': 9})
plt.rcParams.update({'savefig.dpi': 300})
plt.rcParams.update({'savefig.facecolor': 'azure'})
plt.rcParams.update({'savefig.edgecolor': 'k'})
plt.rcParams.update({'savefig.format': 'png'})
plt.rcParams.update({'savefig.jpeg_quality': 95})
plt.rcParams.update({'savefig.pad_inches': 0.05})
plt.rcParams.update({'text.color': 'k'})

# set folder
if _platform == "darwin":
    folder = mac_folder
elif _platform == "win32":
    folder = win_folder
    
# define functions
def mean_angle(deg):
    deg = deg[~deg.isnull()]
    if len(deg) > 2:
        deg = deg.astype(float)
        deg_avg = degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))
        if deg_avg <0:
            deg_avg = 360 + deg_avg
        return deg_avg
    return np.mean(deg)
    
def summarize(df,begin,r_int):
    df_mean = df[df.Time >= begin].resample(r_int, how = 'mean')
    df_mean.rename(columns=lambda x: x+'_avg', inplace=True)
    df_min = df[df.Time >= begin].resample(r_int, how = 'min')
    df_min.rename(columns=lambda x: x+'_min', inplace=True)
    df_max = df[df.Time >= begin].resample(r_int, how = 'max')
    df_max.rename(columns=lambda x: x+'_max', inplace=True)
    df_mean['WindDir'] = df['WindDir'].resample(r_int, how = mean_angle)    
    df = pd.merge(df_mean,pd.merge(df_min,df_max,left_index=True, right_index=True),left_index=True, right_index=True)
    df.drop(['DateUTC_max', 'DateUTC_min', 'Time_max', 'Time_min', 'WindDir_avg',
       'WindDir_max', 'WindDir_min', 'WindDirection_max', 'PrecipDaily_min', 'PrecipHourly_min',
       'WindDirection_min'], axis=1, inplace=True)
    df = df.reindex_axis(sorted(df.columns), axis=1)           
    return df
    
def rawlimit_date(df,begin):
    df = df[df.Time >= begin]
    df = df.reindex_axis(sorted(df.columns), axis=1)           
    return df
    
def rawlimit_daynite(df):
    df_day = df[(df.Time.map(lambda x:x.hour) >= 6) & (df.Time.map(lambda x:x.hour) < 18)]
    df_nit = df[~df.Time.isin(df_day.Time)]
    df_day = df_day.reindex_axis(sorted(df_day.columns), axis=1)
    df_nit = df_nit.reindex_axis(sorted(df_nit.columns), axis=1)           
               
    return df_day,df_nit
    
def rawlimit_season(df):
    df_spr = df[(df.Time.map(lambda x:x.month) >= 3) & (df.Time.map(lambda x:x.month) <= 5)]
    df_smr = df[(df.Time.map(lambda x:x.month) >= 6) & (df.Time.map(lambda x:x.month) <= 8)]
    df_fal = df[(df.Time.map(lambda x:x.month) >= 9) & (df.Time.map(lambda x:x.month) <= 11)]
    df_not = df[(df.Time.map(lambda x:x.month) >= 3) & (df.Time.map(lambda x:x.month) <= 11)]
    df_win = df[~df.Time.isin(df_not.Time)]
    df_win = df_win.reindex_axis(sorted(df_win.columns), axis=1)
    df_spr = df_spr.reindex_axis(sorted(df_spr.columns), axis=1)
    df_smr = df_smr.reindex_axis(sorted(df_smr.columns), axis=1)
    df_fal = df_fal.reindex_axis(sorted(df_fal.columns), axis=1)
          
    return df_win,df_spr,df_smr,df_fal

def treat(df,col1,a,b,col2):
    try:
        df.loc[df[col1] <= a, col1] = b
        df.rename(columns={col1:col2}, inplace=True)
    except:
        pass

def wr(df,p_int,name,latest,city,state,lat,long,elev):
    df.loc[:,'WindDir'] = np.round((df.WindDir / (360 / p_int)),0) * (360 / p_int)
    df.loc[df['WindDir'] == 360,'WindDir'] = 0
    windNA = df[df.Windspeed == 0]
    wind00 = df[(df.Windspeed > 0) & (df.Windspeed <= 1)]
    wind01 = df[(df.Windspeed > 1) & (df.Windspeed <= 2)]
    wind02 = df[(df.Windspeed > 2) & (df.Windspeed <= 5)]
    wind05 = df[(df.Windspeed > 5) & (df.Windspeed <= 10)]
    wind10 = df[df.Windspeed >= 10]

    gb_wind00 = wind00.groupby(['WindDir'])['Windspeed'].count().reset_index().rename(columns={'Windspeed': 'wind00'})
    gb_wind01 = wind01.groupby(['WindDir'])['Windspeed'].count().reset_index().rename(columns={'Windspeed': 'wind01'})
    gb_wind02 = wind02.groupby(['WindDir'])['Windspeed'].count().reset_index().rename(columns={'Windspeed': 'wind02'})
    gb_wind05 = wind05.groupby(['WindDir'])['Windspeed'].count().reset_index().rename(columns={'Windspeed': 'wind05'})
    gb_wind10 = wind10.groupby(['WindDir'])['Windspeed'].count().reset_index().rename(columns={'Windspeed': 'wind10'})
    
    df_deg = pd.DataFrame(data = np.linspace(0, 360, p_int, endpoint=False), columns = ['WindDir'])
    
    df_wr = df_deg.merge(gb_wind00, on=['WindDir'], how='left')\
    .merge(gb_wind01, on=['WindDir'], how='left')\
    .merge(gb_wind02, on=['WindDir'], how='left')\
    .merge(gb_wind05, on=['WindDir'], how='left')\
    .merge(gb_wind10, on=['WindDir'], how='left')
    
    df_wr = df_wr.fillna(0)
    df_wr = df_wr.sort(columns=['WindDir'], axis=0, ascending=True)
    
    total = len(windNA)+len(wind00)+len(wind01)+len(wind02)+len(wind05)+len(wind10)
    
    df_wr['windNA'] = len(windNA)/total*100/p_int
    df_wr.wind00 = df_wr.wind00/total*100
    df_wr.wind01 = df_wr.wind01/total*100
    df_wr.wind02 = df_wr.wind02/total*100
    df_wr.wind05 = df_wr.wind05/total*100
    df_wr.wind10 = df_wr.wind10/total*100
    
    theta = np.linspace(0.0, 2 * np.pi, p_int, endpoint=False)
    width_polar = 2 * np.pi / p_int
    theta = theta - (width_polar/2)
    
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(4.5, 4.5)
    ax = plt.subplot(projection="polar")
    ax.set_axisbelow(True)
    ax.spines['polar'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.bar(theta, df_wr.windNA, width=width_polar, color='#3366ff', edgecolor='none', bottom=0.0)
    ax.bar(theta, df_wr.wind00, width=width_polar, color='#009999', edgecolor='k', bottom=df_wr.windNA)
    ax.bar(theta, df_wr.wind01, width=width_polar, color='#00cc00', edgecolor='k', bottom=df_wr.windNA+df_wr.wind00)
    ax.bar(theta, df_wr.wind02, width=width_polar, color='#bfff00', edgecolor='k', bottom=df_wr.windNA+df_wr.wind00+df_wr.wind01)
    ax.bar(theta, df_wr.wind05, width=width_polar, color='#ffcc00', edgecolor='k', bottom=df_wr.windNA+df_wr.wind00+df_wr.wind01+df_wr.wind02)
    ax.bar(theta, df_wr.wind10, width=width_polar, color='#ffff00', edgecolor='k', bottom=df_wr.windNA+df_wr.wind00+df_wr.wind01+df_wr.wind02+df_wr.wind05)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    fig.text(0.02,0.98,sid+' - '+city+', '+state+'\n'+str(lat)+', '+str(long)+', '+str(elev)+'ft\nWind Rose - '+name,
                         verticalalignment='top',
                         horizontalalignment='left',
                         transform=ax.transAxes)
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax.transAxes)
    fig.savefig(folder+'plots/'+sid+'_wr_'+name+'.'+plt.rcParams['savefig.format'])

def wd(df,p_int,name,size,latest,city,state,lat,long,elev):
    width = 1/24/60*int(df.index.freqstr[:-1])
    df.loc[:,'WindDir'] = np.round((df.WindDir / (360 / p_int)),0) * (360 / p_int)
    df.loc[df['WindDir'] == 360,'WindDir'] = 0
    
    windNA = df[df.Windspeed_avg == 0]
    wind00 = df[(df.Windspeed_avg > 0) & (df.Windspeed_avg <= 1)]
    wind01 = df[(df.Windspeed_avg > 1) & (df.Windspeed_avg <= 2)]
    wind02 = df[(df.Windspeed_avg > 2) & (df.Windspeed_avg <= 5)]
    wind05 = df[(df.Windspeed_avg > 5) & (df.Windspeed_avg <= 10)]
    wind10 = df[df.Windspeed_avg >= 10]
    
    gustNA = df[df.WindspeedGust_max == 0]
    gust00 = df[(df.WindspeedGust_max > 0) & (df.WindspeedGust_max <= 1)]
    gust01 = df[(df.WindspeedGust_max > 1) & (df.WindspeedGust_max <= 2)]
    gust02 = df[(df.WindspeedGust_max > 2) & (df.WindspeedGust_max <= 5)]
    gust05 = df[(df.WindspeedGust_max > 5) & (df.WindspeedGust_max <= 10)]
    gust10 = df[df.WindspeedGust_max >= 10]
    
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(sid+' - '+city+', '+state+': '+str(lat)+', '+str(long)+', '+str(elev)+'ft\nWind Plot - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')
    ax1.bar(gustNA.index, np.zeros(len(gustNA))-height, width, edgecolor='none', color='#3366ff')
    ax1.bar(gust00.index, np.zeros(len(gust00))-height, width, edgecolor='none', color='#009999')
    ax1.bar(gust01.index, np.zeros(len(gust01))-height, width, edgecolor='none', color='#00cc00')
    ax1.bar(gust02.index, np.zeros(len(gust02))-height, width, edgecolor='none', color='#bfff00')
    ax1.bar(gust05.index, np.zeros(len(gust05))-height, width, edgecolor='none', color='#ffcc00')
    ax1.bar(gust10.index, np.zeros(len(gust10))-height, width, edgecolor='none', color='#ffff00')
    ax1.bar(windNA.index, np.zeros(len(windNA))-(height/2), width, edgecolor='none', color='#3366ff')
    ax1.bar(wind00.index, np.zeros(len(wind00))-(height/2), width, edgecolor='none', color='#009999')
    ax1.bar(wind01.index, np.zeros(len(wind01))-(height/2), width, edgecolor='none', color='#00cc00')
    ax1.bar(wind02.index, np.zeros(len(wind02))-(height/2), width, edgecolor='none', color='#bfff00')
    ax1.bar(wind05.index, np.zeros(len(wind05))-(height/2), width, edgecolor='none', color='#ffcc00')
    ax1.bar(wind10.index, np.zeros(len(wind10))-(height/2), width, edgecolor='none', color='#ffff00')
    ax2 = ax1.twinx()
    ax2.plot_date(df.index, df.WindspeedGust_max, marker = '', color='g', linestyle='-', linewidth=1)
    ax2.plot_date(df.index, df.Windspeed_avg, marker = '', color='b', linestyle='-', linewidth=1)
    ylim2max = np.round(((df.WindspeedGust_max.max()+5) / 5),0) * 5
    ylim2min = -ylim2max/4
    ax2.set_ylim([-ylim2min,ylim2max])
    ax2.set_yticks(np.linspace(0,ylim2max,ylim2max/5))
    ax2.set_ylabel('Windspeed / Gust (mph)')
    ax1.plot_date(df.index, df.WindDir, 'rx', ms=size)
    ax1.set_yticks(np.linspace(0,360,(p_int/4)+1))
    ax1.set_ylabel('Wind Direction (deg)')
    ax1.set_ylim([-height,360])
    ax1.grid(b=True, which='both', color='k',linestyle='-')
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax1.transAxes)
    fig.savefig(folder+'plots/'+sid+'_wd_'+name+'.'+plt.rcParams['savefig.format'])

def tdhd(df,name,lw,latest,city,state,lat,long,elev):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(sid+' - '+city+', '+state+': '+str(lat)+', '+str(long)+', '+str(elev)+'ft\nTemp, Dewpoint, Humidity - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')
    ax1.plot_date(df.index, df.Dewpoint_avg, marker = '', color='b', linestyle='-', linewidth=lw)
    ax1.fill_between(df.index, df.Dewpoint_min, df.Dewpoint_max, facecolor='b', alpha=.2)
    ax1.plot_date(df.index, df.Temp_avg, marker = '', color='r', linestyle='-', linewidth=lw)
    ax1.fill_between(df.index, df.Temp_min, df.Temp_max, facecolor='r', alpha=.2)
    ax2 = ax1.twinx()    
    ax2.plot_date(df.index, df.Humidity_avg, marker = '', color='g', linestyle='-', linewidth=lw)
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Temp/Dewpoint (degF)')
    ax2.set_ylabel('Humidity (%)')
    ax1.grid(b=True, which='both', color='k',linestyle='-')
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax1.transAxes)
    fig.savefig(folder+'plots/'+sid+'_tdhd_'+name+'.'+plt.rcParams['savefig.format'])

def pppd(df,name,lw,latest,city,state,lat,long,elev):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(sid+' - '+city+', '+state+': '+str(lat)+', '+str(long)+', '+str(elev)+'ft\nPressure + Precipiation - '+name)
    ax1 = plt.subplot()
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')
    ax1.plot_date(df.index, df.Pressure_avg, marker = '', color='y', linestyle='-', linewidth=lw)
    ax2 = ax1.twinx()
    ax2.plot_date(df.index, df.PrecipHourly_avg, marker = '', color='c', linestyle='-', linewidth=lw)
    ax2.plot_date(df.index, df.PrecipDaily_avg, marker = '', color='b', linestyle='-', linewidth=lw)
    ax1.set_xlabel('Date')
    ax1.set_ylabel('Pressure (inHg)')
    ax2.set_ylabel('Precipitation (in)')
    ax1.grid(b=True, which='both', color='k',linestyle='-')
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax1.transAxes)
    fig.savefig(folder+'plots/'+sid+'_pppd_'+name+'.'+plt.rcParams['savefig.format'])  

def dTdts(df,name,size,latest,city,state,lat,long,elev):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(sid+' - '+city+', '+state+': '+str(lat)+', '+str(long)+', '+str(elev)+'ft\nSolar Radiation, dT/dt, Temp - '+name)
    ax = plt.subplot()
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    axcb = ax.scatter(df.SolarRadiation_avg, df.dTdt_avg, c=df.Temp_avg, vmin=-20, vmax=120, s=size, cmap=plt.cm.jet, alpha=0.75, lw = 0)
    ax.set_xlabel('Solar Radiation (W/m^2)')
    ax.set_ylabel('dT/dt (degF/hr)')
    ax.grid(b=True, which='both', color='k',linestyle='-')
    cbar = plt.colorbar(axcb)
    cbar.set_label('Temp (degF)')
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax.transAxes)
    fig.savefig(folder+'plots/'+sid+'_dTdts_'+name+'.'+plt.rcParams['savefig.format'])
    
def dTdtd(df,name,size,latest,city,state,lat,long,elev):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(sid+' - '+city+', '+state+': '+str(lat)+', '+str(long)+', '+str(elev)+'ft\nTemp + dT/dt - '+name)
    ax = plt.subplot()
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    v = max(abs(df.dTdt_avg.max()),abs(df.dTdt_avg.min()))
    axcb = ax.scatter(df.index, df.Temp_avg, c=df.dTdt_avg, s=size, cmap=plt.cm.jet, vmin=-v, vmax=v, alpha=0.75, lw = 0)
    ax.set_xlim([df.index.min(),df.index.max()])
    ax.set_xlabel('Date')
    ax.set_ylabel('Temp (degF)')
    ax.grid(b=True, which='both', color='k',linestyle='-')
    cbar = plt.colorbar(axcb)
    cbar.set_label('dT/dt (degF/hr)')
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax.transAxes)
    fig.savefig(folder+'plots/'+sid+'_dTdtd_'+name+'.'+plt.rcParams['savefig.format'])    

def tdhs(df,name,size,latest,city,state,lat,long,elev):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 5)
    plt.title(sid+' - '+city+', '+state+': '+str(lat)+', '+str(long)+', '+str(elev)+'ft\nTemp, Dewpoint, Humidity - '+name)
    ax = plt.subplot()
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    axcb = ax.scatter(df.Dewpoint_avg, df.Temp_avg, c=df.Humidity_avg, s=size, cmap=plt.cm.winter, vmin=0, vmax=100, alpha=0.75, lw = 0)
    ax.set_xlabel('Dewpoint (degF)')
    ax.set_ylabel('Temp (degF)')
    ax.grid(b=True, which='both', color='k',linestyle='-')
    cbar = plt.colorbar(axcb)
    cbar.set_label('Humidity (%)')
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax.transAxes)
    fig.savefig(folder+'plots/'+sid+'_tdhs_'+name+'.'+plt.rcParams['savefig.format'])
    
def combo(df,name,lw,latest,city,state,lat,long,elev):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(7.5, 9)
    fig.autofmt_xdate()
    fig.suptitle(sid+' - '+city+', '+state+'\n'+str(lat)+', '+str(long)+', '+str(elev)+'ft\nCombo Plot - '+name)
    ax1 = plt.subplot(5,1,1)
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')
    ax1.plot_date(df.index, df.Dewpoint_avg, marker = '', color='b', linestyle='-', linewidth=lw)
    ax1.fill_between(df.index, df.Dewpoint_min, df.Dewpoint_max, facecolor='b', alpha=.2)
    ax1.plot_date(df.index, df.Temp_avg, marker = '', color='r', linestyle='-', linewidth=lw)
    ax1.fill_between(df.index, df.Temp_min, df.Temp_max, facecolor='r', alpha=.2)
    ax1.set_ylabel('Temp/Dewpoint (degF)')
    ax1.set_ylim([0,120])
    ax1.grid(b=True, which='both', color='k',linestyle='-')
    
    ax2 = plt.subplot(5,1,2)  
    ax2.plot_date(df.index, df.Humidity_avg, marker = '', color='g', linestyle='-', linewidth=lw)
    ax2.fill_between(df.index, df.Humidity_min, df.Humidity_max, facecolor='g', alpha=.2)
    ax2.set_ylabel('Humidity (%)')
    ax22 = ax2.twinx()
    ax22.plot_date(df.index, df.PrecipHourly_max, marker = '', color='c', linestyle='-', linewidth=lw)
    ax22.plot_date(df.index, df.PrecipDaily_max, marker = '', color='b', linestyle='-', linewidth=lw)
    ax22.set_ylabel('Precipitation (in)')
    ax2.set_ylim([0,100])
    ax22.set_ylim([0,2])
    ax2.grid(b=True, which='both', color='k',linestyle='-')
    
    ax3 = plt.subplot(5,1,3)  
    ax3.plot_date(df.index, df.Pressure_avg, marker = '', color='y', linestyle='-', linewidth=lw)
    ax3.fill_between(df.index, df.Pressure_min, df.Pressure_max, facecolor='y', alpha=.2)
    ax3.set_ylabel('Pressure (inHg)')
    ax32 = ax3.twinx()
    ax32.plot_date(df.index, df.CloudBase_avg, marker = '', color='c', linestyle='-', linewidth=lw)
    ax32.fill_between(df.index, df.CloudBase_min, df.CloudBase_max, facecolor='c', alpha=.2)
    ax32.set_ylabel('Minimum Cloudbase (ft)')
    ax3.grid(b=True, which='both', color='k',linestyle='-')
    
    ax4= plt.subplot(5,1,4)
    try:
        ax4.plot_date(df.index, df.SolarRadiation_avg, marker = '', color='orange', linestyle='-', linewidth=lw)
        ax4.fill_between(df.index, df.SolarRadiation_min, df.SolarRadiation_max, facecolor='orange', alpha=.2)    
    except:
        pass
    ax4.set_ylabel('Solar Radiation (W/m^2)')
    ax42 = ax4.twinx()
    ax42.plot_date(df.index, df.dTdt_avg, marker = '', color='r', linestyle='-', linewidth=lw)
    ax42.fill_between(df.index, df.dTdt_min, df.dTdt_max, facecolor='r', alpha=.2)
    ax42.set_ylabel('dT/dt (degF/hr)')
    ax4.grid(b=True, which='both', color='k',linestyle='-')
    
    ax5= plt.subplot(5,1,5)  
    ax5.plot_date(df.index, df.WindDir, 'rx', ms=5*lw)
    ax5.set_yticks(np.linspace(0,360,(p_int/4)+1))
    ax5.set_ylabel('Wind Direction (deg)')
    ax52 = ax5.twinx()
    ax52.plot_date(df.index, df.WindspeedGust_max, marker = '', color='g', linestyle='-', linewidth=lw)
    ax52.plot_date(df.index, df.Windspeed_avg, marker = '', color='b', linestyle='-', linewidth=lw)
    ax52.set_ylabel('Windspeed / Gust (mph)')
    ax5.grid(b=True, which='both', color='k',linestyle='-')
    ax5.set_xlabel('Date')
    
    fig.text(0.98,0.02,latest,fontsize=7,verticalalignment='bottom',horizontalalignment='right',transform=ax5.transAxes)    
    fig.savefig(folder+'plots/'+sid+'_combo_'+name+'.'+plt.rcParams['savefig.format'])

def main(start_date,sid):
    # get station data
    url_sd = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID='+sid
    soup_sd = bs4.BeautifulSoup(urllib.request.urlopen(url_sd))
    #full = soup_sd.find('full').getText()
    neigh = soup_sd.find('neighborhood').getText()
    city = soup_sd.find('city').getText()
    state = soup_sd.find('state').getText()
    lat = float(soup_sd.find('latitude').getText())
    long = float(soup_sd.find('longitude').getText())
    elev = int(soup_sd.find('elevation').getText()[:-3])
    latest = soup_sd.find('observation_time').getText()
    
    # initialize the data pull
    start_date = dt.datetime.strptime(start_date, "%Y-%m-%d").date()
    today = dt.datetime.today().date()
    lim_mon = today+relativedelta(months=-1)
    lim_3mo = today+relativedelta(months=-3)
    lim_wek = today+relativedelta(days=-7)
    fetch_range = pd.DataFrame(columns=['Date'],data=pd.date_range(start_date, today))

    try:
        data = pd.read_csv(folder+'data/'+sid+'_data.csv',index_col=False)
        time.sleep(5)
        data_dates = pd.DataFrame(columns=['Date'],data=pd.to_datetime(data.Time).dt.date.unique())
        recent_range = pd.DataFrame(columns=['Date'],data=pd.date_range(lim_wek, today).date)
        fetch_dates = fetch_range[~fetch_range.Date.isin(data_dates.Date)].append(recent_range).drop_duplicates()
    except:
        fetch_dates = fetch_range

    # print header data
    print('-- WX_SCRAPER --')
    print('(c) Adam Reeder')
    print('station id: '+sid)
    print('neighborhood: '+neigh)
    print('city: '+city)
    print('state: '+state)
    print('lat: '+str(lat))
    print('long: '+str(long))
    print('elevation: '+str(elev))
    print('start date: '+str(start_date))
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
    
        url = 'http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID='+sid+'&day='+str(day)+'&month='+str(month)+'&year='+str(year)+'&graphspan=day&format=1'
        file = folder+'data/'+sid+'_'+year_str+'_'+month_str+'_'+day_str+'_fetch.csv'
        soup = bs4.BeautifulSoup(urllib.request.urlopen(url)).text
        if soup.count('\n') > 1:
            text_file = open(file, 'w')
            text_file.write(soup)
            text_file.close()
            print(year_str+'-'+month_str+'-'+day_str+': fetched')  
        else:
            print(year_str+'-'+month_str+'-'+day_str+': empty')
            continue
    
    # compile to main dataframe
    try:
        for csv1 in glob.glob(folder+'data/'+sid+'*_fetch.csv'):
            data = data.append(pd.read_csv(csv1,index_col=False))
            os.remove(csv1)
    except:
        data = pd.DataFrame()
        for csv1 in glob.glob(folder+'data/'+sid+'*_fetch.csv'):
            data = data.append(pd.read_csv(csv1,index_col=False))
            os.remove(csv1)
            
    data.to_csv(folder+'data/'+sid+'_data.csv',index_col=False)
            
    df = data
    del data
            
    print('data files loaded')
    
    df['Time'] = pd.to_datetime(df.Time)
    df['DateUTC'] = pd.to_datetime(df.DateUTC)
    df.index = df.Time
    df = df.sort(['DateUTC'],ascending = True)
    
    # data treatment
    treat(df,'DewpointF',-50,pd.np.nan,'Dewpoint')
    treat(df,'HourlyPrecipIn',0,0,'PrecipHourly')
    treat(df,'Humidity',0,pd.np.nan,'Humidity')
    treat(df,'PressureIn',0,pd.np.nan,'Pressure')
    treat(df,'SolarRadiationWatts/m^2',-100,pd.np.nan,'SolarRadiation')
    treat(df,'TemperatureF',-100,pd.np.nan,'Temp')
    treat(df,'WindDirectionDegrees',-.01,pd.np.nan,'WindDir')
    treat(df,'WindSpeedGustMPH',-.01,pd.np.nan,'WindspeedGust')
    treat(df,'WindSpeedMPH',-.01,pd.np.nan,'Windspeed')
    treat(df,'dailyrainin',-.01,pd.np.nan,'PrecipDaily')
    
    df.loc[df['Windspeed'] == 0, 'WindDir'] = pd.np.nan
    
    print('data treated')
    
    # calculations
    df['SID'] = sid
    df['Lat'] = lat
    df['Long'] = long
    df['Elev'] = elev
    df['City'] = city
    df['State'] = state
    df['CloudBase'] = (((df['Temp'] - df['Dewpoint']) / 4.4) * 1000) + elev
    df['dT'] = df.Temp - df.Temp.shift(1)
    df['dP'] = df.Pressure - df.Pressure.shift(1)
    df['dt'] = df.Time - df.Time.shift(1)
    df['dt'] = df.dt.astype('timedelta64[s]')/(60*60)
    df['dTdt'] = df['dT'] / df['dt']
    df['dPdt'] = df['dP'] / df['dt']
    df.drop(['dt','dT','dP','SoftwareType'], axis=1, inplace=True)
    
    print('calculations completed')
    
    # archive
    if save_archive == True:
        df.to_csv(folder+'archive/'+sid+'_archive.csv')
        df.drop(['SID','Lat','Long','Elev', 'City', 'State'], axis=1, inplace=True)
        print('data archived')
    else:
        print('data archive not selected')
    
    # create limited dataframes and summaries
    df_rawl_tdy = rawlimit_date(df,today)
    df_rawl_wek = rawlimit_date(df,lim_wek)
    df_rawl_mon = rawlimit_date(df,lim_mon)
    df_rawl_3mo = rawlimit_date(df,lim_3mo)
    df_rawl_all = rawlimit_date(df,start_date)
    
    df_rawl_day,df_rawl_nit = rawlimit_daynite(df)
    df_rawl_win,df_rawl_spr,df_rawl_smr,df_rawl_fal = rawlimit_season(df)
    
    print('timeframes created')
    
    df_summ_tdy = summarize(df,today,str(summ_tdy)+'min')
    df_summ_wek = summarize(df,lim_wek,str(summ_wek)+'min')
    df_summ_mon = summarize(df,lim_mon,str(summ_mon)+'min')
    df_summ_3mo = summarize(df,lim_3mo,str(summ_3mo)+'min')
    df_summ_all = summarize(df,start_date,str(summ_all)+'min')
    
    print('summaries created')

    # create combo plots
    try:
        combo(df_summ_tdy,'day',1,latest,city,state,lat,long,elev)
        combo(df_summ_wek,'week',1,latest,city,state,lat,long,elev)
        combo(df_summ_mon,'month',1,latest,city,state,lat,long,elev)
        #combo(df_summ_3mo,'3mo',1,latest,city,state,lat,long,elev)
        combo(df_summ_all,'all',1,latest,city,state,lat,long,elev)
    except:
        pass
        
    # create windroses
    try:
        wr(df_rawl_tdy,p_int,'day',latest,city,state,lat,long,elev)
        wr(df_rawl_wek,p_int,'week',latest,city,state,lat,long,elev)
        wr(df_rawl_mon,p_int,'month',latest,city,state,lat,long,elev)
        wr(df_rawl_3mo,p_int,'3mo',latest,city,state,lat,long,elev)
        wr(df_rawl_all,p_int,'all',latest,city,state,lat,long,elev)
        wr(df_rawl_day,p_int,'day',latest,city,state,lat,long,elev)
        wr(df_rawl_nit,p_int,'night',latest,city,state,lat,long,elev)
    except:
        pass
    
    # create wind plots
    try:
        wd(df_summ_tdy,p_int,'day',4,latest,city,state,lat,long,elev)
        wd(df_summ_wek,p_int,'week',4,latest,city,state,lat,long,elev)
        wd(df_summ_mon,p_int,'month',4,latest,city,state,lat,long,elev)
        wd(df_summ_3mo,p_int,'3mo',4,latest,city,state,lat,long,elev)
        wd(df_summ_all,p_int,'all',4,latest,city,state,lat,long,elev)
    except:
        pass
    
    # create plots
    try:
        dTdts(df_summ_tdy,'day',75,latest,city,state,lat,long,elev)
        dTdts(df_summ_wek,'week',75,latest,city,state,lat,long,elev)
        dTdts(df_summ_mon,'month',75,latest,city,state,lat,long,elev)
        #dTdts(df_summ_3mo,'3mo',75,latest,city,state,lat,long,elev)
        dTdts(df_summ_all,'all',75,latest,city,state,lat,long,elev)
    except:
        pass
    
    # create plots
    try:
        tdhs(df_summ_tdy,'day',75,latest,city,state,lat,long,elev)
        #tdhs(df_summ_wek,'week',75,latest,city,state,lat,long,elev)
        #tdhs(df_summ_mon,'month',75,latest,city,state,lat,long,elev)
        #tdhs(df_summ_3mo,'3mo',75,latest,city,state,lat,long,elev)
        tdhs(df_summ_all,'all',75,latest,city,state,lat,long,elev)
    except:
        pass
    
    # create plots
    try:
        dTdtd(df_summ_tdy,'day',75,latest,city,state,lat,long,elev)
        dTdtd(df_summ_wek,'week',75,latest,city,state,lat,long,elev)
        #dTdtd(df_summ_mon,'month',75,latest,city,state,lat,long,elev)
        #dTdtd(df_summ_3mo,'3mo',75,latest,city,state,lat,long,elev)
        #dTdtd(df_summ_all,'all',75,latest,city,state,lat,long,elev)
    except:
        pass
        
    print('plots completed')
    
# main loop
while 1<2:
    for sid in sid_list:
        try:
            main(start_date,sid)
        except:
            print(sid+' failed')
    time.sleep(60*sleep_time)    
