# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 14:43:48 2016

@author: areed145
"""

#import matplotlib.pyplot as plt
import pylab as plt
import pandas as pd

plt.figure()
ax=plt.gca()
#ax.axis('off')

size = 8

table = [['431.5', '1049.4', '799.6', '2149.8', '917.9'],
        ['292.2', '717.8', '456.4', '1368.5', '865.6'],
        ['213.8', '636.0', '305.7', '1175.2', '796.0'],
        ['124.6', '555.4', '153.2', '677.2', '192.5'],
        ['66.4', '174.3', '75.1', '577.9', '32.0']]
cols=['Today','Last Day','Last Month','Last 3-Months','Last Year']
rows=['Average Temp','Min Temp','Max Temp','Wind','Test']
table = pd.DataFrame(data=table, columns=cols, index=rows)

cols_num = 1/(len(cols)+1)
rows_num = 1/(len(rows)+1)

for i,col in enumerate(cols):
    for j,row in enumerate(rows):
        plt.text((i+1)*cols_num,j*rows_num,table[col][row],size=size)

for i,col in enumerate(cols):
    plt.text((i+1)*cols_num,(len(rows))*rows_num,cols[i],size=size)
    
for j,row in enumerate(rows):
    plt.text(0,j*rows_num,rows[j],size=size)

plt.show()
plt.savefig("table.png")
