import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# plot time series data for all 8 areas, hourly data for a year
INDIR = sys.argv[1]

file = os.path.join(INDIR, 'Load.CSV')

data = pd.read_csv(file)

# Create a function to convert time t to months
def t_to_month_label(t):
    months = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    ]
    # Days per month in a leap year
    days_per_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    hours_per_month = [days * 24 for days in days_per_month]
    cumulative_hours = np.cumsum(hours_per_month)
    
    for i, hours in enumerate(cumulative_hours):
        if t <= hours:
            return months[i]
    return "Dec"  # Just in case t exceeds 8784

# Create an array of month labels for the x-axis
num_ticks = 12
t_ticks = np.linspace(0, 8784, num_ticks, dtype=int)
month_labels = [t_to_month_label(tick) for tick in t_ticks]

groups = data.groupby('row')

linestyle_tuple = [
    (0, (3, 5, 1, 5)),
    (0, (1, 1)),
    (0, (3, 5, 1, 5)),
    (5, (10, 3)),
    (0, (5, 1)),
    (0, ()),
    (0,(2,2)),
    (0, (5, 10)),
    (0, (3, 10, 1, 10)),
    (0, (3, 1, 1, 1))]

plt.figure(figsize=(12, 6))

# for name, group in groups:
#     plt.plot(group['label'],group['newval'],linestyle=linestyle_tuple[name],label=name)
    
# Calculate and plot the sum and average for all areas
sum_Y = data.groupby('label')['newval'].sum()
avg_Y = data.groupby('label')['newval'].mean()
sum_O = data.groupby('label')['val'].sum()
# print("All loads: ", sum_Y, " AVG: ", sum_Y.mean())
plt.plot(sum_O.index, sum_O.values, label='Orignial Sum', linestyle='-', linewidth=2, color='grey')
plt.plot(sum_Y.index, sum_Y.values, label='Scaled Sum', linestyle=':', linewidth=2, color='deepskyblue')
# plt.plot(avg_Y.index, avg_Y.values, label='Average', linestyle=':', linewidth=2, color='black')

# Customize x-axis labels to display months
plt.xticks(ticks=t_ticks, labels=month_labels, rotation=45,fontsize=20, fontname='Arial')
plt.yticks(fontsize=20, fontname='Arial') 
# plt.title('Load at each area', fontsize=25, fontname='Arial')
plt.ylabel('Load (MW)', fontsize=25, fontname='Arial')
plt.xlabel('Month', fontsize=25, fontname='Arial')
# plt.legend(title='Area')
plt.legend()
plt.show()


