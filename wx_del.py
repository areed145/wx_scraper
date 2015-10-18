# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""

import pandas as pd

obs_comb = pd.read_pickle('/Users/areed145/Dropbox/wx_scraper/obs_data.pk1')

obs_comb = obs_comb.drop(['temp_c','dewpoint_c','pressure_mb','windchill_c','heat_index_c'], axis = 1)

obs_comb['station_id'] = 'KCABAKER38'

obs_comb.to_pickle('/Users/areed145/Dropbox/wx_scraper/obs_data.pk1')