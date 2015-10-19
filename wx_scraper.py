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

sid = 'KCABAKER38'
resample_int = '2min'
width = 1/24/60*2 #interval is a day, need to divide by reample int
p_int = 16
folder = '/Users/areed145/Dropbox/wx_scraper/'

url = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID='+sid
last_obs_time = None

while 2 > 1:
    
    try:
        soup = bs4.BeautifulSoup(urllib.request.urlopen(url))
    except:
        print('WU check failed... waiting')
        time.sleep(15)
        continue

    cur_obs_time = pd.Timestamp(soup.find('observation_time_rfc822').getText()).tz_convert('US/Pacific')
    cur_obs_date = cur_obs_time.strftime('%Y%m%d')
    pickle = folder+'pickles/pickle_'+cur_obs_date+'.pk1'
    archive = folder+'archives/archive_'+cur_obs_date+'.csv'
    plot_obs = folder+'plots/obs_'+cur_obs_date+'.png'

    try:
        df_store = pd.read_pickle(pickle)
        last_obs_time = df_store.index.max()
    except:
        df_store = pd.DataFrame()

    if cur_obs_time != last_obs_time:

        station_id = soup.find('station_id').getText()
        full = soup.find('full').getText()
        neighborhood = soup.find('neighborhood').getText()
        city = soup.find('city').getText()
        state = soup.find('state').getText()
        zipcode = soup.find('zip').getText()
        station_type = soup.find('station_type').getText()
        wind_string = soup.find('wind_string').getText()
        obs_tz = soup.find('observation_time').getText()

        df = pd.DataFrame(index=[cur_obs_time])
        df['station_id'] = soup.find('station_id').getText()
        df['latitude'] = soup.find('latitude').getText()
        df['longitude'] = soup.find('longitude').getText()
        df['elevation'] = soup.find('elevation').getText()
        df['temp_f'] = soup.find('temp_f').getText()
        df['dewpoint_f'] = soup.find('dewpoint_f').getText()
        df['relative_humidity'] = soup.find('relative_humidity').getText()
        df['wind_dir'] = soup.find('wind_dir').getText()
        df['wind_degrees'] = soup.find('wind_degrees').getText()
        df['wind_mph'] = soup.find('wind_mph').getText()
        df['wind_gust_mph'] = soup.find('wind_gust_mph').getText()
        df['pressure_in'] = soup.find('pressure_in').getText()
        df['solar_radiation'] = soup.find('solar_radiation').getText()
        df['uv'] = soup.find('uv').getText()
        df['heat_index_f'] = soup.find('heat_index_f').getText()
        df['windchill_f'] = soup.find('windchill_f').getText()
        df['precip_1hr_in'] = soup.find('precip_1hr_in').getText()
        df['precip_1hr_metric'] = soup.find('precip_1hr_metric').getText()
        df['precip_today_in'] = soup.find('precip_today_in').getText()

        df = df.convert_objects(convert_numeric=True)
        df['cloud_base'] = (((df['temp_f'] - df['dewpoint_f']) / 4.4) * 1000) + 400
        df['obs_d'] = cur_obs_date
        df_store = df_store.append(df)
        df_store.to_pickle(pickle)
        df_store.to_csv(archive)

        last_obs_time = cur_obs_time

        #plotting
        obs_plot_n = df_store.resample(resample_int, how = 'mean')
        obs_plot = df_store.resample(resample_int, how = 'mean').interpolate()        
        obs_plot_max = df_store.resample(resample_int, how = 'max').interpolate()
        
        obs_plot['wind_degrees_u'] = np.round((obs_plot['wind_degrees'] / (360 / p_int)),0) * (360 / p_int)
        obs_plot.loc[obs_plot.wind_degrees_u == 360,'wind_degrees_u'] = 0        
        obs_windNA = obs_plot[obs_plot.wind_mph == 0]
        obs_wind00 = obs_plot[(obs_plot.wind_mph > 0) & (obs_plot.wind_mph <= 2)]
        obs_wind02 = obs_plot[(obs_plot.wind_mph > 2) & (obs_plot.wind_mph <= 5)]
        obs_wind05 = obs_plot[(obs_plot.wind_mph > 5) & (obs_plot.wind_mph <= 10)]
        obs_wind10 = obs_plot[obs_plot.wind_mph >= 10]

        obs_gustNA = obs_plot_max[obs_plot_max.wind_gust_mph == 0]
        obs_gust00 = obs_plot_max[(obs_plot_max.wind_gust_mph > 0) & (obs_plot_max.wind_gust_mph <= 2)]
        obs_gust02 = obs_plot_max[(obs_plot_max.wind_gust_mph > 2) & (obs_plot_max.wind_gust_mph <= 5)]
        obs_gust05 = obs_plot_max[(obs_plot_max.wind_gust_mph > 5) & (obs_plot_max.wind_gust_mph <= 10)]
        obs_gust10 = obs_plot_max[obs_plot_max.wind_gust_mph >= 10]

        plt.close("all")

        fig = mpl.pyplot.gcf()
        fig.set_size_inches(10, 25)
        plt.suptitle(station_id+': '+full+'\n'+cur_obs_date+', Last Observation: '+last_obs_time.strftime('%H:%M:%S'), fontsize=20)
        
        #plt.subplot(7, 1, 1)
        plt.subplot2grid((10,2), (0, 0), colspan=2)
        plt.bar(obs_plot.index, obs_plot['cloud_base']/100, width, color='#99CCFF', edgecolor='none')
        plt.ylabel('Min Cloudbase (ft MSL)')
        #plt.ylim([0,100])
        plt.grid(b=True, which='both', color='k',linestyle='-')

        #plt.subplot(7, 1, 2)
        plt.subplot2grid((10,2), (1, 0), colspan=2)
        plt.plot(obs_plot.index, obs_plot['temp_f'], color='#3366FF',linestyle='-')
        plt.plot(obs_plot.index, obs_plot['dewpoint_f'], color='#CCFF33',linestyle='-')

        try:
            plt.plot(obs_plot_n.index, obs_plot_n['heat_index_f'], color='#FF0066',linestyle='-')
            plt.plot(obs_plot_n.index, obs_plot_n['wind_chill_f'], color='#00FFCC',linestyle='-')
        except:
            pass

        plt.ylabel('Temp/Dewpt (F)')
        #plt.ylim([0,110])
        plt.grid(b=True, which='both', color='k',linestyle='-')

        #plt.subplot(7, 1, 3)
        plt.subplot2grid((10,2), (2, 0), colspan=2)
        plt.plot(obs_plot.index, obs_plot['pressure_in'], color='#FF5C33',linestyle='-')
        plt.ylabel('Pres (in)')
        #plt.ylim([29.6,30.6])
        plt.grid(b=True, which='both', color='k',linestyle='-')

        #plt.subplot(7, 1, 4)
        plt.subplot2grid((10,2), (3, 0), colspan=2)
        plt.bar(obs_gustNA.index, obs_gustNA['wind_degrees']-obs_gustNA['wind_degrees']+359, width, edgecolor='none', color='#FFFFFF')
        plt.bar(obs_gust00.index, obs_gust00['wind_degrees']-obs_gust00['wind_degrees']+359, width, edgecolor='none', color='#FFCC00')
        plt.bar(obs_gust02.index, obs_gust02['wind_degrees']-obs_gust02['wind_degrees']+359, width, edgecolor='none', color='#FF9900')
        plt.bar(obs_gust05.index, obs_gust05['wind_degrees']-obs_gust05['wind_degrees']+359, width, edgecolor='none', color='#FF0000')
        plt.bar(obs_gust10.index, obs_gust10['wind_degrees']-obs_gust10['wind_degrees']+359, width, edgecolor='none', color='#FF3385')
        plt.bar(obs_windNA.index, obs_windNA['wind_degrees']-obs_windNA['wind_degrees']+250, width, edgecolor='none', color='#FFFFFF')
        plt.bar(obs_wind00.index, obs_wind00['wind_degrees']-obs_wind00['wind_degrees']+250, width, edgecolor='none', color='#FFCC00')
        plt.bar(obs_wind02.index, obs_wind02['wind_degrees']-obs_wind02['wind_degrees']+250, width, edgecolor='none', color='#FF9900')
        plt.bar(obs_wind05.index, obs_wind05['wind_degrees']-obs_wind05['wind_degrees']+250, width, edgecolor='none', color='#FF0000')
        plt.bar(obs_wind10.index, obs_wind10['wind_degrees']-obs_wind10['wind_degrees']+250, width, edgecolor='none', color='#FF3385')
        plt.plot(obs_plot.index, obs_plot['wind_degrees'], 'k-', ms=2)
        plt.ylabel('Wind Dir (deg)')
        plt.ylim([0,359])
        plt.grid(b=True, which='both', color='k',linestyle='-')        

        #plt.subplot(7, 1, 5)
        plt.subplot2grid((10,2), (4, 0), colspan=2)
        plt.bar(obs_plot.index, obs_plot['solar_radiation']/100, width, color='#FFDB4D', edgecolor='none')
        plt.bar(obs_plot.index, obs_plot['uv'], width, color='#DB4DB8', edgecolor='none')
        plt.ylabel('Solar Rad (w/mm^2) / UV')
        plt.ylim([0,10])
        plt.grid(b=True, which='both', color='k',linestyle='-')

        #plt.subplot(7, 1, 6)
        plt.subplot2grid((10,2), (5, 0), colspan=2)
        plt.bar(obs_plot.index, obs_plot['relative_humidity'], width, color='#256EB8', edgecolor='none')
        plt.bar(obs_plot.index, 100 - obs_plot['relative_humidity'], width, bottom=obs_plot['relative_humidity'], color='#B2CCFF', edgecolor='none')
        plt.ylabel('Rel Humidity (%)')
        plt.ylim([0,100])
        plt.grid(b=True, which='both', color='k',linestyle='-')

        #plt.subplot(7, 1, 7)
        plt.subplot2grid((10,2), (6, 0), colspan=2)
        plt.bar(obs_plot.index, obs_plot['precip_today_in'], width, color='#335CD6', edgecolor='none')
        plt.plot(obs_plot.index, obs_plot['precip_1hr_in'], color='#00CCFF',linestyle='-')
        plt.ylabel('Precip (in)')
        plt.ylim([0,0.5])
        plt.grid(b=True, which='both', color='k',linestyle='-')

        windNA = obs_windNA.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'windNA'})
        wind00 = obs_wind00.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind00'})
        wind02 = obs_wind02.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind02'})
        wind05 = obs_wind05.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind05'})
        wind10 = obs_wind10.groupby(['wind_degrees_u'])['wind_degrees'].count().reset_index().rename(columns={'wind_degrees': 'wind10'})
        
        df = windNA.merge(wind00, on=['wind_degrees_u'], how='outer')\
                    .merge(wind02, on=['wind_degrees_u'], how='outer')\
                    .merge(wind05, on=['wind_degrees_u'], how='outer')\
                    .merge(wind10, on=['wind_degrees_u'], how='outer')
                    
        df = df.fillna(0)
        df = df.sort(columns=['wind_degrees_u'], axis=0, ascending=True)
                    
        total = df.windNA.sum()+df.wind00.sum()+df.wind02.sum()+df.wind05.sum()+df.wind10.sum()
        
        df['calm'] = df.windNA.sum()/total*100/16
        df.wind00 = df.wind00/total*100
        df.wind02 = df.wind02/total*100
        df.wind05 = df.wind05/total*100
        df.wind10 = df.wind10/total*100
                    
        theta = np.linspace(0.0, 2 * np.pi, p_int, endpoint=False)
        width_polar = 2 * np.pi / p_int
        theta = theta - (width/2)
        
        ax = plt.subplot2grid((10,2), (7, 0), colspan=2, rowspan=3, projection="polar")
        #ax = plt.subplot(111, polar=True)
        bars = ax.bar(theta, df.calm, width=width_polar, color='#8AB8E6', edgecolor='none', bottom=0.0)
        bars = ax.bar(theta, df.wind00, width=width_polar, color='#FFCC00', edgecolor='#FFFFFF', bottom=df.calm)
        bars = ax.bar(theta, df.wind02, width=width_polar, color='#FF9900', edgecolor='#FFFFFF', bottom=df.calm+df.wind00)
        bars = ax.bar(theta, df.wind05, width=width_polar, color='#FF0000', edgecolor='#FFFFFF', bottom=df.calm+df.wind00+df.wind02)
        bars = ax.bar(theta, df.wind10, width=width_polar, color='#FF3385', edgecolor='#FFFFFF', bottom=df.calm+df.wind00+df.wind02+df.wind05)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)

        fig.savefig(plot_obs, dpi=200)
        
        print('WU checked: '+str(cur_obs_time))

    time.sleep(15)