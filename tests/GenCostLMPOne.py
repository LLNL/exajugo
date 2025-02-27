import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from collections import defaultdict
from scipy import stats

INDIR = sys.argv[1]
Day = [1,35]
BRatio = [0.0,0.1,0.3]
# Day = [1,35,80,93,129,176,195,355]
# BRatio = [0.0,0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3]

# Dictionary to store data organized by BRatio values
dataGCost = {}
all_L_values = defaultdict(lambda: defaultdict(list))
average_L_over_time = {}

# Function to extend the list for a specific row r and dimension d
def extend_list(r, d):
    all_L_values[r][d].extend(df['LMP'].tolist())
    
# Populate all_L_values
for r in BRatio:
    r_dict = {}
    for d in Day:
        r_dict[d] = []
    all_L_values.append(r_dict)
    
# Create list of file names and read data 
for r in BRatio:
    dfsG = []
    l_sum_over_t = None
    count_d = 0
    
    r_key = (r,)
    for d in Day:
        file = os.path.join(INDIR,'BusSol_24h_'+'NOZ'+str(d)+'D'+str(r)+'.CSV')
        
        if os.path.exists(file):
            df = pd.read_csv(file, index_col=0)  # Use first column as index (periods)
            dfsG.extend(df['GencostAvg'].values)
    
        # Append L values to the corresponding r key in all_L_values dictionary
            if d not in all_L_values[r]:
                extend_list(r, d)
            #     all_L_values[r][d] = []
            # all_L_values[r][d].extend(df['LMP'].tolist())
            
         # Sum L values over periods t for the current (r, d) pair
            if l_sum_over_t is None:
                l_sum_over_t = df.groupby('time')['LMP'].sum()
            else:
                l_sum_over_t += df.groupby('time')['LMP'].sum()

            count_d += 1

    if count_d > 0:
        avg_l_over_t = l_sum_over_t / count_d
        average_L_over_time[r] = avg_l_over_t
        
    if dfsG:
        dataGCost[r_key] = np.array(dfsG)
                
                
# print(dfsG)
# Plot histograms for LMP values for each r
plt.figure(figsize=(14, 10))
num_rows = int(np.ceil(len(BRatio) / 2))
for i, r in enumerate(all_L_values):
    all_L_values_r = []
    for d in all_L_values[r]:
        all_L_values_r.extend(all_L_values[r][d])
    plt.subplot(num_rows, 2, i + 1)
    plt.hist(all_L_values_r, bins=20, alpha=0.7, label=f'R = {r}', color='blue')
    plt.xlabel('LMP ($/MWh)', fontname='Arial')
    plt.ylabel('Frequency', fontname='Arial')
    plt.legend()
    plt.grid(True)
    # plt.xlim(min(all_L_values[r]), 25) 

plt.tight_layout()
# plt.show()

# Plotting average L vs time for each r
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

plt.figure(figsize=(10, 6))

i = 0
for r in average_L_over_time:
    avg_l_over_t = average_L_over_time[r]
    plt.plot(range(24), avg_l_over_t.values, linestyle=linestyle_tuple[i], linewidth=3, marker='o', label=f'R = {r}')
    i += 1

plt.xlabel('Time (hour)', fontsize=25, fontname='Arial')
plt.ylabel('Average LMPs ($/MWh)', fontsize=25, fontname='Arial')
plt.title('Average LMPs at different bidding ratio, single-bus setting', fontsize=25, fontname='Arial')
plt.xticks(fontsize=20, fontname='Arial') 
plt.yticks(fontsize=20, fontname='Arial')
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.show()
        
# Calculate the mean, standard deviation, and range for each BRatio for each day
# Initialize dictionaries to store statistics
means = defaultdict(lambda: defaultdict(float))
stds = defaultdict(lambda: defaultdict(float))
ranges = defaultdict(lambda: defaultdict(float))
conf_intervals = defaultdict(lambda: defaultdict(tuple))

overall_means = {}
overall_stds = {}
overall_ranges = {}
overall_conf_intervals = {}

for r in BRatio:
    for d in Day:
        n = len(all_L_values[r][d])
        means[r][d] = np.mean(all_L_values[r][d])
        stds[r][d] = np.std(all_L_values[r][d],ddof=1)
        conf_intervals[r][d] = stats.t.interval(0.95, n-1, loc=means[r][d], scale=stds[r][d]/np.sqrt(n))
        ranges[r][d] = (np.min(all_L_values[r][d]), np.max(all_L_values[r][d]))
    
    n = len(all_L_values_r)
    overall_means[r] = np.mean(all_L_values[r])
    overall_stds[r] = np.std(all_L_values[r])
    overall_conf_intervals[r] = stats.t.interval(0.95, n-1, loc=overall_means[r], scale=overall_stds[r]/np.sqrt(n))
    overall_ranges[r] = (np.min(all_L_values[r]), np.max(all_L_values[r]))
   
print("ranges for each r and d:", ranges)
print("overall ranges:", overall_ranges)       
# Plotting
plt.figure(figsize=(14, 10))
colors = plt.cm.tab20(np.linspace(0, 1, len(Day)))
# Plot mean values with error bars and range bands for each d
for d_index, d in enumerate(Day):
    mean_vals = [means[r][d] for r in BRatio]
    std_vals = [stds[r][d] for r in BRatio]
    range_vals = [ranges[r][d] for r in BRatio]
    lower_conf_vals = [conf_intervals[r][d][0] for r in BRatio]
    upper_conf_vals = [conf_intervals[r][d][1] for r in BRatio]

    plt.plot(BRatio, mean_vals, label=f'Mean for {d}',color=colors[d_index],linestyle=linestyle_tuple[d_index],linewidth=3)
    # print(mean_vals)
    plt.errorbar(BRatio, mean_vals, yerr=[np.array(mean_vals) - np.array(lower_conf_vals), np.array(upper_conf_vals) - np.array(mean_vals)], fmt='o',color=colors[d_index])
    plt.fill_between(BRatio, [r[0] for r in range_vals], [r[1] for r in range_vals],color=colors[d_index], alpha=0.1)

# Plot overall mean values with error bars and range bands
overall_mean_vals = [overall_means[r] for r in BRatio]
overall_std_vals = [overall_stds[r] for r in BRatio]
overall_range_vals = [overall_ranges[r] for r in BRatio]
overall_lower_conf_vals = [overall_conf_intervals[r][0] for r in BRatio]
overall_upper_conf_vals = [overall_conf_intervals[r][1] for r in BRatio]

plt.plot(BRatio, overall_mean_vals, label='Overall Mean', color='red', linestyle='-')
plt.errorbar(BRatio, overall_mean_vals,yerr=[np.array(overall_mean_vals) - np.array(overall_lower_conf_vals), np.array(overall_upper_conf_vals) - np.array(overall_mean_vals)], fmt='o', color='red', capsize=5)
plt.fill_between(BRatio, [r[0] for r in overall_range_vals], [r[1] for r in overall_range_vals], color='red', alpha=0.1)

# aggregated_data = {}
# for r_key, g_values in dataGCost.items():
#     r_value = r_key[0]  # Extract float value from tuple
    
#     mean_g = np.mean(g_values)
#     std_g = np.std(g_values)
#     min_g = np.min(g_values)
#     max_g = np.max(g_values)
    
#     aggregated_data[r_key] = {
#         'mean': mean_g,
#         'std': std_g,
#         'min': min_g,
#         'max': max_g,
#         'all_g_values': g_values 
#     }

# # plot mean lines, standard deviation bands, and range bands
# plt.figure(figsize=(14, 10))

# for r_key, stats in aggregated_data.items():
#     r_value = r_key[0]  # Extract float value from tuple
    
#     mean_g = stats['mean']
#     std_g = stats['std']
#     min_g = stats['min']
#     max_g = stats['max']
#     all_g_values = stats['all_g_values']
    
#     # Plot data line with points 
#     plt.plot([r_value], [mean_g], marker='o', markersize=8, color='blue', linewidth=5)
#     # Error bars 
#     plt.errorbar([r_value], [mean_g], yerr=std_g, fmt='', color='blue', capsize=5)
#     # Fill between (blue shaded region with transparency)
#     plt.fill_between([r_value], min_g, max_g, color='blue', alpha=0.3)
#      # Scatter plot of all G values (gray)
#     plt.scatter([r_value] * len(all_g_values), all_g_values, color='gray', alpha=0.5)  # Scatter plot of all G values

# Customize plot
plt.xlabel('Bidding ratio', fontsize=25, fontname='Arial')
plt.ylabel('Average gen cost ($/MWh)', fontsize=25, fontname='Arial')
plt.title('Average generation cost at differnt bidding ratio, single-bus setting', fontsize=25, fontname='Arial')
plt.legend()
plt.xticks(BRatio)
# plt.ylim(16.65,17.15)
plt.xticks(fontsize=20, fontname='Arial') 
plt.yticks(fontsize=20, fontname='Arial') 
plt.xticks(rotation=60)  # Rotate x-axis labels for better readability
# Format x-axis tick labels to display float values with minimal precision
plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.3g}'))
plt.tight_layout()
plt.show()



                
        
    


