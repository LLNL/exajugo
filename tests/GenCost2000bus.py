import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from collections import defaultdict
from scipy import stats

# new format
INDIR = sys.argv[1]
# Day = [224]
# BRatio = [0.0,0.01,0.025,0.05,0.075,0.1,0.13,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]

Day = [1,35,87,95,110,129,176,195,224,268,290,300,330,335,355]
BRatio = [0.0,0.002,0.01,0.025,0.05,0.075,0.1,0.13]

# Initialize all_L_values as a list of dictionaries
all_L_values = defaultdict(lambda: defaultdict(list))
    
# Extend the lists with actual data
for r in BRatio:
    for d in Day:
        # for generation cost
        file = os.path.join(INDIR,'Loss_24h_'+'M10OZ'+str(d)+'D'+str(r)+'.xlsx')
        
        if os.path.exists(file):
            df = pd.read_excel(file, sheet_name='Loss', nrows=24, index_col=0)  # Use first column as index (periods)
            df['GencostAvg'] = df['GenCost']/df['GenTot']
            # Append L values to the corresponding r key in all_L_values dictionary
            if d not in all_L_values[r]:
                all_L_values[r][d] = []
            all_L_values[r][d].extend(df['GencostAvg'].tolist())  
        else:
            print("file doesn't exist", file)
        # # ## for LMPs
        # file = os.path.join(INDIR,'GLPAllbus_24h_'+'TNOZ'+str(d)+'D'+str(r)+'.CSV')  
          
        # if os.path.exists(file):
        #     df = pd.read_csv(file)
        #     # Append L values to the corresponding r key in all_L_values dictionary
        #     if d not in all_L_values[r]:
        #         all_L_values[r][d] = []
        #     all_L_values[r][d].extend(df['LMP'].tolist())  
        # else:
        #     print("file doesn't exist", file)    
     
# print(len(all_L_values[0.0][1]))         
# Initialize dictionaries to store statistics
means = defaultdict(lambda: defaultdict(float))
ranges = defaultdict(lambda: defaultdict(float))
conf_intervals = defaultdict(lambda: defaultdict(tuple))

# Calculate statistics for each r and d
for r in BRatio:
    means[r] = {}
    conf_intervals[r] = {}
    ranges[r] = {}
    for d in Day:
        values = all_L_values[r][d]
        n = len(values)
        mean = np.mean(values)
        std = np.std(values, ddof=1)
        conf_interval = stats.t.interval(0.95, n-1, loc=mean, scale=std/np.sqrt(n))
        means[r][d] = mean
        conf_intervals[r][d] = conf_interval
        ranges[r][d] = (np.min(values), np.max(values))
        
# Calculate overall statistics for each r across all d
overall_means = {}
overall_conf_intervals = {}
overall_ranges = {}

for r in BRatio:
    all_values_for_r = []
    for d in Day:
        all_values_for_r.extend(all_L_values[r][d])
    n = len(all_values_for_r)
    mean = np.mean(all_values_for_r)
    std = np.std(all_values_for_r, ddof=1)
    conf_interval = stats.t.interval(0.95, n-1, loc=mean, scale=std/np.sqrt(n))
    overall_means[r] = mean
    overall_conf_intervals[r] = conf_interval
    overall_ranges[r] = (np.min(all_values_for_r), np.max(all_values_for_r))

# Plotting
plt.figure(figsize=(12, 8))

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
    (5, (10, 3)),
    (0, (5, 1)),
    (0, ()),
    (0,(2,2)),
    (0, (5, 10)),
    (0, (3, 10, 1, 10)),
    (0, (3, 1, 1, 1))]

# Generate a list of unique colors
colors = plt.cm.tab20(np.linspace(0, 1, len(Day)))

# Plot mean values with 95% confidence intervals and range bands for each d
for d_index, d in enumerate(Day):
    mean_vals = [means[r][d] for r in BRatio]
    # print(mean_vals)
    lower_conf_vals = [conf_intervals[r][d][0] for r in BRatio]
    upper_conf_vals = [conf_intervals[r][d][1] for r in BRatio]
    range_vals = [ranges[r][d] for r in BRatio]

    # plt.plot(BRatio, mean_vals, label=f'Day {d}', color=colors[0],linestyle=linestyle_tuple[1],linewidth=3)
    # plt.errorbar(BRatio, mean_vals, yerr=[np.array(mean_vals) - np.array(lower_conf_vals), np.array(upper_conf_vals) - np.array(mean_vals)], fmt='o', color=colors[0], capsize=5)
    # plt.fill_between(BRatio, [r[0] for r in range_vals], [r[1] for r in range_vals], color=colors[0], alpha=0.1)

# Plot overall mean values with 95% confidence intervals and range bands
overall_mean_vals = [overall_means[r] for r in BRatio]
overall_lower_conf_vals = [overall_conf_intervals[r][0] for r in BRatio]
overall_upper_conf_vals = [overall_conf_intervals[r][1] for r in BRatio]
overall_range_vals = [overall_ranges[r] for r in BRatio]

print(overall_mean_vals)

plt.plot(BRatio, overall_mean_vals, label='Average', color=colors[0],linestyle=linestyle_tuple[1],linewidth=5)
plt.errorbar(BRatio, overall_mean_vals, yerr=[np.array(overall_mean_vals) - np.array(overall_lower_conf_vals), np.array(overall_upper_conf_vals) - np.array(overall_mean_vals)], fmt='o', color=colors[0], capsize=5)
# plt.fill_between(BRatio, [r[0] for r in overall_range_vals], [r[1] for r in overall_range_vals], color='red', alpha=0.1)
plt.fill_between(BRatio, [r[0] for r in overall_range_vals], [r[1] for r in overall_range_vals], color=colors[0], alpha=0.08)

# Customize plot
plt.xlabel('bidding ratio', fontsize=34, fontname='Arial')
plt.ylabel('average gen cost ($/MWh)', fontsize=34, fontname='Arial')
# plt.title('Average generation cost at different bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')
# plt.ylabel('Average LMPs ($/MWh)', fontsize=25, fontname='Arial')
# plt.title('Average LMPs at different bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')
# plt.legend(fontsize=15, ncol=4)
plt.xticks(BRatio)
# Manually set uniform grid lines using plt.gca() (get current axis)
# plt.gca().yaxis.set_major_locator(plt.MultipleLocator(0.1))  # Uniform grid spacing of 0.1 for y-axis
plt.gca().xaxis.set_major_locator(plt.MultipleLocator(0.01))    # Uniform grid spacing of 1 for x-axis
plt.grid(True,linestyle='--')

plt.ylim(16,26)
# plt.ylim(16.6,17.2)
plt.xticks(fontsize=35, fontname='Arial') 
plt.yticks(fontsize=35, fontname='Arial') 
plt.xticks(rotation=60)  # Rotate x-axis labels for better readability
# Format x-axis tick labels to display float values with minimal precision
plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.3g}'))
plt.tight_layout()
plt.show()