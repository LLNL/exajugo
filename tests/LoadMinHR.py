import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# plot time series data for all 8 areas, hourly data for a year
INDIR = sys.argv[1]

file = os.path.join(INDIR, 'Load.CSV')

data = pd.read_csv(file)

# Group by hour and sum the N values for each hour
hourly_data = data.groupby('label')['newval'].sum().reset_index()
hourly_data['day'] = ((hourly_data['label'] -1)// 24) + 1
# print(hourly_data.head())
# Find the five hours with the minimum total N load
top_5_min_hours = hourly_data.nsmallest(5, 'newval')
top_5_max_hours = hourly_data.nlargest(5, 'newval')

# Print the results
print(f"The first five hours with the minimum total are:\n{top_5_min_hours}")
print(f"The first five hours with the max total are:\n{top_5_max_hours}")

hourly_data_sort = hourly_data.sort_values(by='newval')
# outfile = os.path.join(INDIR, 'HourlyLoad.xlsx')
# hourly_data_sort.to_excel(outfile)

# hourly_data_fea = hourly_data_sort[hourly_data_sort['newval'] >= 24564.95]
# hourly_data_fea_sort = hourly_data_fea.sort_values(by=['day','label'])
# outfile_fea = os.path.join(INDIR, 'HourlyLoadFea.xlsx')
# hourly_data_fea_sort.to_excel(outfile_fea)

hourly_data_infea = hourly_data_sort[hourly_data_sort['newval'] < 24564.95]
hourly_data_infea_sort = hourly_data_infea.sort_values(by=['day','label'])
outfile_infea = os.path.join(INDIR, 'HourlyLoadinFea.xlsx')
hourly_data_infea_sort.to_excel(outfile_infea)

# # Print the results
# print("The five hours with the minimum total load are:")
# for index, row in top_5_min_hours.iterrows():
#     print(f"Hour {row['label']}: Total load = {row['newval']}")

# print("The five hours with the maximum total load are:")
# for index, row in top_5_max_hours.iterrows():
#     print(f"Hour {row['label']}: Total load = {row['newval']}")


# # Display the data for these hours
# print("\nThe data for these hours is:")
# print(top_5_min_hours)




