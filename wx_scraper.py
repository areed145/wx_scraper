# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""
import bs4
import urllib
import pandas as pd
import threading

obs_time = None
obs_single = pd.DataFrame()
obs_comb = pd.DataFrame()

def check_wu():
  threading.Timer(5.0, check_wu).start()
  url = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID=KCABAKER38'
  soup = bs4.BeautifulSoup(urllib.request.urlopen(url))  
  print ('wu_checked')  
  return soup

#soup = check_wu()

url = 'http://api.wunderground.com/weatherstation/WXCurrentObXML.asp?ID=KCABAKER38'
soup = bs4.BeautifulSoup(urllib.request.urlopen(url)) 

observation_time_rfc822 = soup.find('observation_time_rfc822').getText()

if observation_time_rfc822 != obs_time:
    obs_single['observation_time_rfc822'] = soup.find('observation_time_rfc822').getText()
    obs_single['full'] = soup.find('full').getText()
    obs_single['neighborhood'] = soup.find('neighborhood').getText()
    obs_single['city'] = soup.find('city').getText()
    obs_single['state'] = soup.find('state').getText()
    obs_single['zip'] = soup.find('zip').getText()
    obs_single['latitude'] = soup.find('latitude').getText()
    obs_single['longitude'] = soup.find('longitude').getText()
    obs_single['elevation'] = soup.find('elevation').getText()
    obs_single['station_id'] = soup.find('station_id').getText()
    obs_single['station_type'] = soup.find('station_type').getText()
    obs_single['temp_f'] = soup.find('temp_f').getText()
    obs_single['temp_c'] = soup.find('temp_c').getText()
    obs_single['relative_humidity'] = soup.find('relative_humidity').getText()
    obs_single['wind_string'] = soup.find('wind_string').getText()
    obs_single['wind_dir'] = soup.find('wind_dir').getText()
    obs_single['wind_degrees'] = soup.find('wind_degrees').getText()
    obs_single['wind_mph'] = soup.find('wind_mph').getText()
    obs_single['wind_gust_mph'] = soup.find('wind_gust_mph').getText()
    obs_single['pressure_mb'] = soup.find('pressure_mb').getText()
    obs_single['pressure_in'] = soup.find('pressure_in').getText()
    obs_single['dewpoint_f'] = soup.find('dewpoint_f').getText()
    obs_single['dewpoint_c'] = soup.find('dewpoint_c').getText()
    obs_single['heat_index_f'] = soup.find('heat_index_f').getText()
    obs_single['heat_index_c'] = soup.find('heat_index_c').getText()
    obs_single['windchill_f'] = soup.find('windchill_f').getText()
    obs_single['windchill_c'] = soup.find('windchill_c').getText()
    obs_single['solar_radiation'] = soup.find('solar_radiation').getText()
    obs_single['uv'] = soup.find('uv').getText()
    obs_single['precip_1hr_in'] = soup.find('precip_1hr_in').getText()
    obs_single['precip_1hr_metric'] = soup.find('precip_1hr_metric').getText()
    obs_single['precip_today_in'] = soup.find('precip_today_in').getText()
    
    if len(obs_comb) < 1:
        obs_comb = obs_single
    else:
        obs_comb = obs_comb.append(obs_single)
        
    obs_time = soup.find('observation_time_rfc822').getText()
    
  



