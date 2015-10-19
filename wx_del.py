# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:20:13 2015

@author: areed145
"""

import pandas as pd
import numpy as np

p_int = 22.5

df = pd.read_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/pickle_20151018.pk1')

df['wind_degrees_u'] = np.round((df['wind_degrees'] / p_int),0)*p_int
df.loc[df.wind_degrees_u == 360,'wind_degrees_u'] = 0

df.to_pickle('/Users/areed145/Dropbox/wx_scraper/pickles/pickle_20151018.pk1')