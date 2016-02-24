# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import bs4
import urllib
import pandas as pd
import datetime as dt
import glob
#import matplotlib as mpl
import matplotlib.pyplot as plt
from cmath import rect, phase
from math import radians, degrees
import numpy as np
from dateutil.relativedelta import relativedelta
#from matplotlib.dates import DateFormatter
#import pytz

#####################################################################################################

# input data
start_date = '2015-01-01'
sid = 'KCABAKER38'
#folder = 'C:/Users/bvjs/Python/python-3.4.3.amd64/bvjs/wx_scraper/'
folder = '/Users/areed145/Dropbox/wx_scraper/'
p_int = 32
width = 1/24/60*5
height = 90

# initialize the data pull
start_date = dt.datetime.strptime(start_date, "%Y-%m-%d").date()
today = dt.datetime.today().strftime("%Y-%m-%d")
today = dt.datetime.strptime(today, "%Y-%m-%d").date()
M = today+relativedelta(months=-1)
TM = today+relativedelta(months=-3)
W = today+relativedelta(days=-7)
delta = (today - start_date).days

# station elevation
url_elev = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID='+sid
soup_elev = bs4.BeautifulSoup(urllib.request.urlopen(url_elev))
elev = soup_elev.find('elevation').getText()
elev = int(elev[:-3])

#####################################################################################################

# define functions
def mean_angle(deg):
    if len(deg) > 2:
        deg = deg.astype(float)
        deg_avg = degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))
        if deg_avg <0:
            deg_avg = 360 + deg_avg
        return deg_avg
    return np.mean(deg)
    
def limiter(df,begin,r_int):
    df_mean = df[df.Time >= begin].resample(r_int, how = 'mean')
    df_mean.rename(columns=lambda x: x+'_avg', inplace=True)
    df_min = df[df.Time >= begin].resample(r_int, how = 'min')
    df_min.rename(columns=lambda x: x+'_min', inplace=True)
    df_max = df[df.Time >= begin].resample(r_int, how = 'max')
    df_max.rename(columns=lambda x: x+'_max', inplace=True)
    df_mean['WindDir'] = df['WindDir'].resample(r_int, how = mean_angle)
    
    df = pd.merge(df_mean,pd.merge(df_min,df_max,left_index=True, right_index=True),left_index=True, right_index=True)
    df.drop(['Clouds_avg', 'Clouds_max', 'Clouds_min', 'Conditions_avg', 'Conditions_max',
       'Conditions_min', 'DateUTC_max', 'DateUTC_min', 'Elevation_avg', 'Elevation_max',
       'Elevation_min', 'SoftwareType_max', 'SoftwareType_min', 'Time_max', 'Time_min',
       'WindDir_avg', 'WindDir_max', 'WindDir_min', 'WindDirection_max', 'PrecipDaily_min',
       'PrecipHourly_min', 'WindDirection_min'], axis=1, inplace=True)
    df = df.reindex_axis(sorted(df.columns), axis=1)   
        
    return df

def treat(col1,a,b,col2):
    try:
        df.loc[df[col1] <= a, col1] = b
        df.rename(columns={col1:col2}, inplace=True)
    except:
        pass

def wind(df,p_int,width,name,size):
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

    gb_windNA = windNA.groupby(['WindDir'])['Windspeed_avg'].count().reset_index().rename(columns={'Windspeed_avg': 'windNA'})
    gb_wind00 = wind00.groupby(['WindDir'])['Windspeed_avg'].count().reset_index().rename(columns={'Windspeed_avg': 'wind00'})
    gb_wind01 = wind01.groupby(['WindDir'])['Windspeed_avg'].count().reset_index().rename(columns={'Windspeed_avg': 'wind01'})
    gb_wind02 = wind02.groupby(['WindDir'])['Windspeed_avg'].count().reset_index().rename(columns={'Windspeed_avg': 'wind02'})
    gb_wind05 = wind05.groupby(['WindDir'])['Windspeed_avg'].count().reset_index().rename(columns={'Windspeed_avg': 'wind05'})
    gb_wind10 = wind10.groupby(['WindDir'])['Windspeed_avg'].count().reset_index().rename(columns={'Windspeed_avg': 'wind10'})
    
    df_deg = pd.DataFrame(data = np.linspace(0, 360, p_int, endpoint=False), columns = ['WindDir'])
    
    df_wr = df_deg.merge(gb_windNA, on=['WindDir'], how='left')\
    .merge(gb_wind00, on=['WindDir'], how='left')\
    .merge(gb_wind01, on=['WindDir'], how='left')\
    .merge(gb_wind02, on=['WindDir'], how='left')\
    .merge(gb_wind05, on=['WindDir'], how='left')\
    .merge(gb_wind10, on=['WindDir'], how='left')
    
    df_wr = df_wr.fillna(0)
    df_wr = df_wr.sort(columns=['WindDir'], axis=0, ascending=True)
    
    total = df_wr.windNA.sum()+df_wr.wind00.sum()+df_wr.wind01.sum()+df_wr.wind02.sum()+df_wr.wind05.sum()+df_wr.wind10.sum()
    
    df_wr.windNA = df_wr.windNA.sum()/total*100/p_int
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
    fig.set_size_inches(10, 10)
    ax = plt.subplot(projection="polar")
    ax.set_axisbelow(True)
    ax.spines['polar'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.bar(theta, df_wr.windNA, width=width_polar, color='#3366ff', edgecolor='none', bottom=0.0)
    ax.bar(theta, df_wr.wind00, width=width_polar, color='#009999', edgecolor='none', bottom=df_wr.windNA)
    ax.bar(theta, df_wr.wind01, width=width_polar, color='#00cc00', edgecolor='none', bottom=df_wr.windNA+df_wr.wind00)
    ax.bar(theta, df_wr.wind02, width=width_polar, color='#bfff00', edgecolor='none', bottom=df_wr.windNA+df_wr.wind00+df_wr.wind01)
    ax.bar(theta, df_wr.wind05, width=width_polar, color='#ffcc00', edgecolor='none', bottom=df_wr.windNA+df_wr.wind00+df_wr.wind01+df_wr.wind02)
    ax.bar(theta, df_wr.wind10, width=width_polar, color='#ffff00', edgecolor='none', bottom=df_wr.windNA+df_wr.wind00+df_wr.wind01+df_wr.wind02+df_wr.wind05)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    fig.savefig(folder+'plots/'+sid+'_wr_'+name+'.png', dpi=400)
    
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(15, 5)
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
    ax2.set_ylim([ylim2min,ylim2max])
    ax2.set_yticks(np.linspace(0,ylim2max,(p_int/4)+1))
    ax1.plot_date(df.index, df.WindDir, 'rx', ms=size)
    ax1.set_yticks(np.linspace(0,360,(p_int/4)+1))
    ax1.set_ylabel('Wind Dir (deg)')
    ax1.set_ylim([-height,360])
    ax1.grid(b=True, which='both', color='k',linestyle='-')
    fig.savefig(folder+'plots/'+sid+'_wd_'+name+'.png', dpi=400)

def tdhd(df,name,lw):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(15, 5)
    ax = plt.subplot()
    #ax.suptitle('Last Observation: ')
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.fill_between(df.index, df.Dewpoint_min, df.Dewpoint_max, facecolor='b', alpha=.1)
    ax.plot_date(df.index, df.Dewpoint_avg, marker = '', color='b', linestyle='-', linewidth=lw)
    ax.fill_between(df.index, df.Temp_min, df.Temp_max, facecolor='r', alpha=.1)
    ax.plot_date(df.index, df.Temp_avg, marker = '', color='r', linestyle='-', linewidth=lw)
    #ax.plot_date(df.index, df.Humidity_avg, marker = '', color='g', linestyle='-', linewidth=lw)
    #plt.ylabel('dT/dt (degF/hr)')
    ax.grid(b=True, which='both', color='k',linestyle='-')
    fig.savefig(folder+'plots/'+sid+'_tdhd_'+name+'.png', dpi=400)

def pppd(df,name,lw):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(15, 5)
    ax1 = plt.subplot()
    #plt.suptitle('Last Observation: ')
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(True)
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')
    ax1.fill_between(df.index, df.Pressure_min, df.Pressure_max, facecolor='y', alpha=.1)
    ax1.plot_date(df.index, df.Pressure_avg, marker = '', color='y', linestyle='-', linewidth=lw)
    ax2 = ax1.twinx()
    ax2.plot_date(df.index, df.PrecipHourly_avg, marker = '', color='c', linestyle='-', linewidth=lw)
    ax2.plot_date(df.index, df.PrecipDaily_avg, marker = '', color='b', linestyle='-', linewidth=lw)
    #plt.ylabel('dT/dt (degF/hr)')
    ax1.grid(b=True, which='both', color='k',linestyle='-')
    fig.savefig(folder+'plots/'+sid+'_pppd_'+name+'.png', dpi=400)  

def dTdts(df,name,size):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    ax = plt.subplot()
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.scatter(df.SolarRadiation_avg, df.Temp_avg, c=df.dTdt, s=size, cmap=plt.cm.jet, alpha=0.75, lw = 0)
    #plt.ylabel('dT/dt (degF/hr)')
    ax.grid(b=True, which='both', color='k',linestyle='-')
    fig.savefig(folder+'plots/'+sid+'_dTdts_'+name+'.png', dpi=400)
    
def dTdtd(df,name,size):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(15, 5)
    ax = plt.subplot()
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.scatter(df.index, df.Temp_avg, c=df.dTdt_avg, s=size, cmap=plt.cm.jet, alpha=0.75, lw = 0)
    ax.scatter(df.index, df.Temp_min, c=df.dTdt_min, s=size, cmap=plt.cm.jet, alpha=0.75, lw = 0)
    ax.scatter(df.index, df.Temp_max, c=df.dTdt_max, s=size, cmap=plt.cm.jet, alpha=0.75, lw = 0)
    ax.set_xlim([df.index.min(),df.index.max()])
    #plt.ylabel('dT/dt (degF/hr)')
    ax.grid(b=True, which='both', color='k',linestyle='-')
    fig.savefig(folder+'plots/'+sid+'_dTdtd_'+name+'.png', dpi=400)    

def tdhs(df,name,size):
    #formatter = DateFormatter('%H:%M', tz=pytz.timezone('US/Pacific'))
    plt.close("all")
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    ax = plt.subplot()
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.scatter(df.Dewpoint_avg, df.Temp_avg, c=df.Humidity, s=size, cmap=plt.cm.winter, alpha=0.75, lw = 0)
    #plt.ylabel('dT/dt (degF/hr)')
    ax.grid(b=True, which='both', color='k',linestyle='-')
    fig.savefig(folder+'plots/'+sid+'_tdhs_'+name+'.png', dpi=400)  

#####################################################################################################

# data pull loop
for date_count in range(delta+1):

    date = start_date+dt.timedelta(days=date_count)
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
    file = folder+'data/'+sid+'_'+year_str+'_'+month_str+'_'+day_str+'.csv'

    if date_count == delta:
        try:
            soup = bs4.BeautifulSoup(urllib.request.urlopen(url)).text  
            text_file = open(file, 'w')
            text_file.write(soup)
            text_file.close()  
        except:
            pass
    else:
        try:
            test = pd.read_csv(file,index_col=False)
            del test
        except:
            try:
                soup = bs4.BeautifulSoup(urllib.request.urlopen(url)).text  
                text_file = open(file, 'w')
                text_file.write(soup)
                text_file.close()  
            except:
                pass

#####################################################################################################

# compile to main dataframe  
df = pd.DataFrame()
for csv1 in glob.glob(folder+'data/'+sid+'*.csv'):
    df = df.append(pd.read_csv(csv1,index_col=False))

df['Time'] = pd.to_datetime(df.Time)
df['DateUTC'] = pd.to_datetime(df.DateUTC)
df.index = df.Time
df = df.sort(['DateUTC'],ascending = True)

# data treatment
treat('DewpointF',-50,pd.np.nan,'Dewpoint')
treat('HourlyPrecipIn',0,0,'PrecipHourly')
treat('Humidity',0,pd.np.nan,'Humidity')
treat('PressureIn',0,pd.np.nan,'Pressure')
treat('SolarRadiationWatts/m^2',-100,pd.np.nan,'SolarRadiation')
treat('TemperatureF',-100,pd.np.nan,'Temp')
treat('WindDirectionDegrees',-.01,pd.np.nan,'WindDir')
treat('WindSpeedGustMPH',-.01,pd.np.nan,'WindspeedGust')
treat('WindSpeedMPH',-.01,pd.np.nan,'Windspeed')
treat('dailyrainin',-.01,pd.np.nan,'PrecipDaily')

# calculations
df['Elevation'] = elev
df['CloudBase'] = (((df['Temp'] - df['Dewpoint']) / 4.4) * 1000) + elev
df['dT'] = df.Temp - df.Temp.shift(1)
df['dP'] = df.Pressure - df.Pressure.shift(1)
df['dt'] = df.Time - df.Time.shift(1)
df['dt'] = df.dt.astype('timedelta64[s]')/(60*60)
df['dTdt'] = df['dT'] / df['dt']
df['dPdt'] = df['dP'] / df['dt']
df.drop(['dt','dT','dP'], axis=1, inplace=True)

# archive
df.to_csv(folder+'archive/'+sid+'.csv')

# create limited dataframes
df_today = limiter(df,today,'15min')
df_week = limiter(df,W,'30min')
df_month = limiter(df,M,'120min')
df_3month = limiter(df,TM,'720min')
df_all = limiter(df,start_date,'1440min')

#####################################################################################################

# create windroses
try:
    wind(df_today,p_int,1/24/60*15,'day',4)
    wind(df_week,p_int,1/24/60*30,'week',4)
    wind(df_month,p_int,1/24*2,'month',4)
    wind(df_3month,p_int,1/24*12,'3month',4)
    wind(df_all,p_int,1,'all',4)
except:
    pass

# create date plots
try:
    tdhd(df_today,'day',3)
    tdhd(df_week,'week',2)
    tdhd(df_month,'month',1)
    tdhd(df_3month,'3month',1)
    tdhd(df_all,'all',1)
except:
    pass

# create date plots
try:
    pppd(df_today,'day',3)
    pppd(df_week,'week',2)
    pppd(df_month,'month',1)
    pppd(df_3month,'3month',1)
except:
    pass

# create date plots
try:
    dTdts(df_today,'day',100)
    dTdts(df_week,'week',100)
    dTdts(df_month,'month',100)
except:
    pass

## create date plots
#try:
#    tdhs(df_today,'day',100)
#    tdhs(df_week,'week',100)
#    tdhs(df_month,'month',100)
#except:
#    pass

# create date plots
try:
    dTdtd(df_today,'day',100)
    dTdtd(df_week,'week',75)
    dTdtd(df_month,'month',50)
except:
    pass