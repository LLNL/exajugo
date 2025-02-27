import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# plot time series data for all 8 areas, hourly data for a year
INDIR = sys.argv[1]

file = os.path.join(INDIR, 'Load.CSV')

data = pd.read_csv(file)
data['day'] = ((data['label'] -1)// 24) + 1

# print("1st day: ", min(data['day']), " last day: ", max(data['day']))
# Group by 'day' and sum the values of 'newval' for each day
daily_data = data.groupby('day')['newval'].sum().reset_index()
# print(daily_data.head())
# Sort the days by the total N and get the first five days with the minimum total N
top_5_min_days = daily_data.nsmallest(5, 'newval')
top_5_max_days = daily_data.nlargest(20, 'newval')

# Print the results
print(f"The first five days with the minimum total are:\n{top_5_min_days}")
print(f"The first five days with the max total are:\n{top_5_max_days}")

daily_data_sorted = daily_data.sort_values(by='newval')
daily_data_sorted = daily_data_sorted.reset_index(drop=True)
daily_data_sorted.insert(0,'rank',daily_data_sorted.index)
# daily_data_feasible = daily_data_sorted[daily_data_sorted['rank'] > 90]
# daily_data_feasible_sort = daily_data_feasible.sort_values(by='day')

# outfile = os.path.join(INDIR, 'DailyLoadFea.xlsx')
# daily_data_feasible_sort.to_excel(outfile)

daily_data_infeasible = daily_data_sorted[daily_data_sorted['rank'] <= 90]
daily_data_infeasible_sort = daily_data_infeasible.sort_values(by='day')

outfile_in = os.path.join(INDIR, 'DailyLoadinFea.xlsx')
daily_data_infeasible_sort.to_excel(outfile_in)


# print(daily_data_sorted.iloc[80:100])
# print(daily_data_sorted.head(50))
# Filter and print the data for the 24 hours of each of these days
# for day in top_5_min_days['day']:
#     min_day_data = data[data['day'] == day]
#     print(f"\nThe data N for the 24 hours of day {day} is:")
    # print(min_day_data[['label', 'newval']])

# ## filter to find days that have loads smaller than total min generation capacity
# gen_min = 24*24564.95 
# daily_data_small = daily_data[daily_data['newval'] < gen_min]
# print("Daily loads smaller than generation lower limits:")
# print(daily_data_small['day'])
