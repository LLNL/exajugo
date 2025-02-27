import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from collections import defaultdict
from scipy import stats

# new format
INDIR = sys.argv[1]
# Day = [1]
# Day = [1,35,87,95,129,176,195,224,268,290,330,355]
# Day = [1,35,62,87,95,129,176,195,224,268,290,330,355]
# Day = [1,35]
Day = [1,35,87,95,110,129,176,195,224,268,290,300,330,335,355]
BRatio = [0.002,0.01,0.05,0.1,0.13]
# BRatio = [0.0,0.1,0.3]
# Day = [1,176,195,223,268,355]
# BRatio = [0.0,0.002,0.01,0.025,0.05,0.075,0.1]

# Day = [1,35,164,176,185,195,223,268,355]
# BRatio = [0.0,0.001,0.01,0.025,0.05,0.1]
# BRatio = [0.0,0.002,0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3]

# Dictionary to store data organized by BRatio values
all_L_values = {}
average_L_over_time = {}

# Create list of file names and read data 
for r in BRatio:
    dfsL = []
    l_sum_over_t = None
    count_d = 0
    
    for d in Day:
        file = os.path.join(INDIR,'BusSol_24h_'+'M10OZ'+str(d)+'D'+str(r)+'.CSV')
        
        if os.path.exists(file):
            df = pd.read_csv(file) 
            # Append G values to the corresponding r key in all_L_values dictionary
            if r not in all_L_values:
                all_L_values[r] = []
            all_L_values[r].extend(df['plantpower'].tolist())
            
            # Sum L values over periods t for the current (r, d) pair
            if l_sum_over_t is None:
                l_sum_over_t = df.groupby('time')['plantpower'].sum()
            else:
                l_sum_over_t += df.groupby('time')['plantpower'].sum()

            count_d += 1
        else:
            print("file doesn't exist", file)
        
    if count_d > 0:
        avg_l_over_t = l_sum_over_t / count_d
        average_L_over_time[r] = avg_l_over_t
        
     
# print(len(all_L_values[0.0]))   
# Plot histograms for each r value
# plt.figure(figsize=(6,10))
# num_rows = int(np.ceil(len(all_L_values)))
# for i, r in enumerate(all_L_values):
#     plt.subplot(num_rows, 1, i + 1)
#     mean_value = np.mean(all_L_values[r])
#     # print(mean_value)
#     plt.hist(all_L_values[r], bins=40, alpha=0.7, label=f'R = {r}', color='blue')
#     plt.axvline(mean_value, color='r',linestyle='-',linewidth=3)
#     plt.text(mean_value, -0.05, f'{mean_value:.2f}', va='top', ha='center', color='r')
#     plt.xlabel('LMP ($/MWh)', fontname='Arial')
#     plt.ylabel('Frequency', fontname='Arial')
#     plt.legend()
#     plt.grid(True)
#     # plt.xlim(0, 22)
#     # plt.ylim(0, 100) 
#     # plt.xlim(0,18.5)
#     # plt.ylim(0,200)
#     # plt.xlim(15.5,18.5)

# plt.figtext(0.5, -0.0, f'Days {Day}', ha='center', va='bottom', fontsize=10, color='gray')
# plt.tight_layout()
# # plt.show()

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

plt.figure(figsize=(12, 8))

i = 0
for r in average_L_over_time:
    avg_l_over_t = average_L_over_time[r]
    plt.plot(range(1,25), avg_l_over_t.values, linestyle=linestyle_tuple[i], linewidth=5, marker='o', label=f'$R_b$ = {r}')
    i += 1

plt.xlabel('time (hour)', fontsize=35, fontname='Arial')
plt.ylabel('plant power (MW)', fontsize=35, fontname='Arial')
# plt.title('Average bidding plant power at different bidding ratio, single-bus setting', fontsize=35, fontname='Arial')
plt.xticks(np.arange(0,25,4),fontsize=35, fontname='Arial') 

plt.tick_params(axis='both', which='major', labelsize=35)
plt.yticks(fontsize=35, fontname='Arial')
plt.xticks(fontsize=35, fontname='Arial')
# plt.ylim(13,21)
# plt.ylim(8,16)
plt.legend(fontsize=25)
plt.grid(True)
# plt.figtext(0.5, 0.01, f'Days {Day}', ha='center', va='bottom', fontsize=20, color='gray')
plt.tight_layout()
plt.show()

