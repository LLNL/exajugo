import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# plot time series data for all 8 areas, hourly data for a year
INDIR = sys.argv[1]

file = os.path.join(INDIR, 'case.CSV')

data = pd.read_csv(file)

data['day'] = ((data['label'] -1)// 24) + 1
data = data.rename(columns={'newval': 'val'})

hourly_data = data.groupby('label').agg({'val': 'sum', 'day':'first'}).reset_index()
# print(hourly_data.head())

# # val_filter = 25000
# val_filter = 24565
# val_scale = 1.15
# dayScale = hourly_data.loc[hourly_data['val'] < val_filter, 'day'].unique()

# data['newval'] = data['val']
# data.loc[data['day'].isin(dayScale),'newval'] = data['val'] * val_scale

all_scale = 1.1
data['newval'] = data['val'] * all_scale


outfile = os.path.join(INDIR, 'Load.CSV')
data.to_csv(outfile)