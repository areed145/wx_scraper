# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import bs4
import urllib
import pandas as pd
import time
import matplotlib as mpl
import matplotlib.pyplot as plt

last_obs_time = None

while 2 > 1:
    
    url = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID=KCABAKER38'
    soup = bs4.BeautifulSoup(urllib.request.urlopen(url))

    cur_obs_time = pd.Timestamp(soup.find('observation_time_rfc822').getText()).tz_convert('US/Pacific')
    cur_obs_date = cur_obs_time.strftime('%Y%m%d')
    pickle = '/Users/areed145/Dropbox/wx_scraper/pickles/pickle_'+cur_obs_date+'.pk1'
    archive = '/Users/areed145/Dropbox/wx_scraper/archives/archive_'+cur_obs_date+'.csv'
    plot_obs = '/Users/areed145/Dropbox/wx_scraper/plots/obs_'+cur_obs_date+'.png'
    
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

    obs_plot = df_store
    
    obs_plot = obs_plot.resample('1min', how = 'mean')
    
    plt.close("all")
    
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(10.5, 16.5)
    
    plt.subplot(9, 1, 1)
    plt.title(station_id+': '+full)
    plt.plot(obs_plot.index, obs_plot['temp_f'], 'b-')
    plt.plot(obs_plot.index, obs_plot['dewpoint_f'], 'g-')
    
    try:
        plt.plot(obs_plot.index, obs_plot['heat_index_f'], 'r-')
        plt.plot(obs_plot.index, obs_plot['wind_chill_f'], 'c-')
    except:
        pass
    
    plt.ylabel('Temp/Dewpt (F)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    
    plt.subplot(9, 1, 2)
    plt.plot(obs_plot.index, obs_plot['pressure_in'], 'k-')
    plt.ylabel('Pres (in)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    
    plt.subplot(9, 1, 3)
    plt.plot(obs_plot.index, obs_plot['wind_mph'], 'y-')
    plt.plot(obs_plot.index, obs_plot['wind_gust_mph'], 'r.')
    plt.ylabel('Wind/Gust (mph)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    
    obs_windNA = obs_plot[obs_plot.wind_mph == 0]
    obs_wind00 = obs_plot[(obs_plot.wind_mph > 0) & (obs_plot.wind_mph <= 2)]
    obs_wind02 = obs_plot[(obs_plot.wind_mph > 2) & (obs_plot.wind_mph <= 5)]
    obs_wind05 = obs_plot[obs_plot.wind_mph >= 5]   
    
    plt.subplot(9, 1, 4)
    plt.plot(obs_windNA.index, obs_windNA['wind_degrees'], '.', ms=3, color='#9191B5')
    plt.plot(obs_wind00.index, obs_wind00['wind_degrees'], '.', ms=5, color='#FFCC00')
    plt.plot(obs_wind02.index, obs_wind02['wind_degrees'], '.', ms=6, color='#FF9900')
    plt.plot(obs_wind05.index, obs_wind05['wind_degrees'], '.', ms=6, color='#FF0000')
    plt.ylabel('Wind Dir (deg)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    plt.ylim([0,359])
    
    plt.subplot(9, 1, 5)
    plt.plot(obs_plot.index, obs_plot['cloud_base'], 'c-')
    plt.ylabel('Min Cloudbase (ft MSL)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    
    plt.subplot(9, 1, 6)
    plt.plot(obs_plot.index, obs_plot['uv'], 'm-')
    plt.ylabel('UV')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    
    plt.subplot(9, 1, 7)
    plt.plot(obs_plot.index, obs_plot['solar_radiation'])
    plt.ylabel('Solar Rad (w/m^2)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    
    plt.subplot(9, 1, 8)
    plt.plot(obs_plot.index, obs_plot['relative_humidity'], 'r-')
    #plt.stackplot(obs_plot.index, obs_plot['relative_humidity'], obs_plot['solar_radiation'])
    plt.ylabel('Rel Humidity (%)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')
    
    plt.subplot(9, 1, 9)
    plt.plot(obs_plot.index, obs_plot['precip_today_in'], 'b-')
    plt.ylabel('Precip (in)')
    plt.grid(b=True, which='both', color='0.65',linestyle='-')

    fig.savefig(plot_obs, dpi=200)
    
    print('WU checked: '+str(cur_obs_time))
    
    time.sleep(15)