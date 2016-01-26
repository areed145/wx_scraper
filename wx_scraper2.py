# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import bs4
import csv
import urllib
import pandas as pd
import datetime as dt
import glob
import matplotlib.pyplot as plt

start_date = '2016-01-24'
start_date = dt.datetime.strptime(start_date, "%Y-%m-%d").date()

today = dt.datetime.today().strftime("%Y-%m-%d")
today = dt.datetime.strptime(today, "%Y-%m-%d").date()

delta = (today - start_date).days

sid = 'KCABAKER38'
folder = '/Users/areed145/Documents/GitHub/wx_scraper/data/'

for date_count in range(delta+1):
    
    date = start_date+dt.timedelta(days=date_count)
    year = date.year
    month = date.month
    day = date.day
    
    url = 'http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID='+sid+'&day='+str(day)+'&month='+str(month)+'&year='+str(year)+'&graphspan=day&format=1'
    
    try:
        soup = bs4.BeautifulSoup(urllib.request.urlopen(url)).text
        file = folder+sid+str(year)+str(month)+str(day)+'.csv'
        text_file = open(file, 'w')
        text_file.write(soup)
        text_file.close()
        a = pd.read_csv(file,index_col=False)
    except:
        pass
    
df = pd.DataFrame()

for csv1 in glob.glob(folder+'*.csv'):
    df = df.append(pd.read_csv(csv1,index_col=False))
    
plt.plot(df['DewpointF'].values)
plt.plot(df['TemperatureF'].values)
plt.plot(df['Humidity'].values)
plt.ylim(0,100)