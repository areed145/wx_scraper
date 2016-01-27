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
import matplotlib as mpl
import matplotlib.pyplot as plt
from cmath import rect, phase
from math import radians, degrees
import numpy as np
from dateutil.relativedelta import relativedelta
#from matplotlib.dates import DateFormatter

# input data
start_date = '2015-01-01'
sid = 'KCABAKER38'
folder = 'C:/Users/bvjs/Python/python-3.4.3.amd64/bvjs/wx_scraper/'
folder = '/Users/areed145/Documents/GitHub/wx_scraper/'
p_int = 32

# define functions
def mean_angle(deg):
            if len(deg) > 2:
                deg = deg.astype(float)
                deg_avg = degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))
                if deg_avg <0:
                    deg_avg = 360 + deg_avg
                return deg_avg
            return pd.np.nan
            
def treat(col1,a,b,col2):
    try:
        df.loc[df[col1] <= a, col1] = b
        df.rename(columns={col1:col2}, inplace=True)
    except:
        pass
    
def windrose(df,p_int,name):
    df['WindDir_U'] = np.round((df['WindDir'] / (360 / p_int)),0) * (360 / p_int)
    df.loc[df['WindDir_U'] == 360,'WindDir_U'] = 0
    windNA = df[df.Windspeed == 0].groupby(['WindDir_U'])['WindDir'].count().reset_index().rename(columns={'WindDir': 'windNA'})
    wind00 = df[(df.Windspeed > 0) & (df.Windspeed <= 1)].groupby(['WindDir_U'])['WindDir'].count().reset_index().rename(columns={'WindDir': 'wind00'})
    wind01 = df[(df.Windspeed > 1) & (df.Windspeed <= 2)].groupby(['WindDir_U'])['WindDir'].count().reset_index().rename(columns={'WindDir': 'wind01'})
    wind02 = df[(df.Windspeed > 2) & (df.Windspeed <= 5)].groupby(['WindDir_U'])['WindDir'].count().reset_index().rename(columns={'WindDir': 'wind02'})
    wind05 = df[(df.Windspeed > 5) & (df.Windspeed <= 10)].groupby(['WindDir_U'])['WindDir'].count().reset_index().rename(columns={'WindDir': 'wind05'})
    wind10 = df[df.Windspeed >= 10].groupby(['WindDir_U'])['WindDir'].count().reset_index().rename(columns={'WindDir': 'wind10'})

    df_deg = pd.DataFrame(data = np.linspace(0, 360, p_int, endpoint=False), columns = ['WindDir_U'])
    
    df_wr = df_deg.merge(windNA, on=['WindDir_U'], how='left')\
                .merge(wind00, on=['WindDir_U'], how='left')\
                .merge(wind01, on=['WindDir_U'], how='left')\
                .merge(wind02, on=['WindDir_U'], how='left')\
                .merge(wind05, on=['WindDir_U'], how='left')\
                .merge(wind10, on=['WindDir_U'], how='left')
    
    df_wr = df_wr.fillna(0)
    df_wr = df_wr.sort(columns=['WindDir_U'], axis=0, ascending=True)
    
    total = df_wr.windNA.sum()+df_wr.wind00.sum()+df_wr.wind01.sum()+df_wr.wind02.sum()+df_wr.wind05.sum()+df_wr.wind10.sum()
    
    df_wr['calm'] = df_wr.windNA.sum()/total*100/p_int
    df_wr.windNA = df_wr.windNA/total*100
    df_wr.wind00 = df_wr.wind00/total*100
    df_wr.wind01 = df_wr.wind01/total*100
    df_wr.wind02 = df_wr.wind02/total*100
    df_wr.wind05 = df_wr.wind05/total*100
    df_wr.wind10 = df_wr.wind10/total*100
    
    theta = np.linspace(0.0, 2 * np.pi, p_int, endpoint=False)
    width_polar = 2 * np.pi / p_int
    theta = theta - (width_polar/2)
        
    plt.close("all")
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(10, 10)
    ax = plt.subplot(projection="polar")
    ax.set_axisbelow(True)
    ax.spines['polar'].set_visible(True)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.bar(theta, df_wr.calm, width=width_polar, color='#3366ff', edgecolor='none', bottom=0.0)
    ax.bar(theta, df_wr.wind00, width=width_polar, color='#009999', edgecolor='k', bottom=df_wr.calm)
    ax.bar(theta, df_wr.wind01, width=width_polar, color='#00cc00', edgecolor='k', bottom=df_wr.calm+df_wr.wind00)
    ax.bar(theta, df_wr.wind02, width=width_polar, color='#bfff00', edgecolor='k', bottom=df_wr.calm+df_wr.wind00+df_wr.wind01)
    ax.bar(theta, df_wr.wind05, width=width_polar, color='#ffcc00', edgecolor='k', bottom=df_wr.calm+df_wr.wind00+df_wr.wind01+df_wr.wind02)
    ax.bar(theta, df_wr.wind10, width=width_polar, color='#ffff00', edgecolor='k', bottom=df_wr.calm+df_wr.wind00+df_wr.wind01+df_wr.wind02+df_wr.wind05)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    fig.savefig(folder+'plots/'+sid+'_WR_'+name+'.png', dpi=400)

# initialize the data pull
start_date = dt.datetime.strptime(start_date, "%Y-%m-%d").date()
today = dt.datetime.today().strftime("%Y-%m-%d")
today = dt.datetime.strptime(today, "%Y-%m-%d").date()
M = today+relativedelta(months=-1)
W = today+dt.timedelta(days=-7)
delta = (today - start_date).days

# station elevation
url_elev = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID='+sid
soup_elev = bs4.BeautifulSoup(urllib.request.urlopen(url_elev))
elev = soup_elev.find('elevation').getText()
elev = int(elev[:-3])

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

# combine to main dataframe    
df = pd.DataFrame()
for csv1 in glob.glob(folder+'data/'+sid+'*.csv'):
    df = df.append(pd.read_csv(csv1,index_col=False))

df['Time'] = pd.to_datetime(df.Time)
df['DateUTC'] = pd.to_datetime(df.DateUTC)
df = df.sort(['DateUTC'],ascending = True)
    
# data treatment
treat('DewpointF',-50,pd.np.nan,'Dewpoint')
treat('HourlyPrecipIn',0,0,'PrecipHourly')
treat('Humidity',0,pd.np.nan,'Humidity')
treat('PressureIn',0,pd.np.nan,'Pressure')
treat('SolarRadiation',-100,pd.np.nan,'SolarRadiation')
treat('TemperatureF',-100,pd.np.nan,'Temp')
treat('WindDirectionDegrees',-.01,pd.np.nan,'WindDir')
treat('WindSpeedGustMPH',-.01,pd.np.nan,'WindspeedGust')
treat('WindSpeedMPH',-.01,pd.np.nan,'Windspeed')
treat('dailyrainin',-.01,pd.np.nan,'PrecipDaily')

# calculations
df['Elevation'] = elev
df['CloudBase'] = (((df['Temp'] - df['Dewpoint']) / 4.4) * 1000) + elev
df['dT'] = df.Temp - df.Temp.shift(1)
df['dt'] = df.Time - df.Time.shift(1)
df['dt'] = df.dt.astype('timedelta64[s]')/(60*60)
df['dT/dt'] = df['dT'] / df['dt']
df.drop(['dt','dT'], axis=1, inplace=True)

# archive
df.to_csv(folder+'archive/'+sid+'.csv')

# create windroses
windrose(df[df.Time >= today],p_int,'D')
windrose(df[df.Time >= W],p_int,'W')
windrose(df[df.Time >= M],p_int,'M')
windrose(df,p_int,'A')