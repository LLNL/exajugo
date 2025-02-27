import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# new format
INDIR = sys.argv[1]
# Day = [1,35] 
# Day = [1,35,87,95,129,176,195,224,268,290,330,355]
# Day = [1,35,62,87,95,129,176,195,224,268,290,330,355]
Day = [1,35,87,95,110,129,176,195,268,290,300,330,335,355]
# Day = [1,35,87,95,110,129,176,195,224,268,290,300,330,335,355]
BRatio = np.array([0.0,0.002,0.01,0.05,0.1,0.13])
# Day = [1,35,93,120,130,176,195,224,268,290,305,355]
# BRatio = np.array([0.0,0.002,0.01,0.025,0.05,0.075,0.1,0.14])
# Day = [1,35,164,176,185,195,223,268,355]
# BRatio = np.array([0.0,0.001,0.01,0.025,0.05,0.1,0.15,0.2])

# Create list of file names for two sets (LMPs and GenCost)
file_names_G = {}
# file_names_C = {}
for r in BRatio:
    for d in Day:
        file_names_G[(r,d)] = os.path.join(INDIR,'GLPAllbus_24h_'+'M10OZ'+str(d)+'D'+str(r)+'.CSV')
        # file_names_C[(r,d)] = os.path.join(INDIR,'Loss_24h_'+'NZ'+str(d)+'D'+str(r)+'.xlsx')


# Initialize a list to store G values for each r
all_G_values = {}

# Iterate through each (r, d) pair
for (r, d), file_name in file_names_G.items():
    if os.path.exists(file_name):
        df = pd.read_csv(file_name)
    
        # Append G values to the corresponding r key in all_G_values dictionary
        if r not in all_G_values:
            all_G_values[r] = []
        all_G_values[r].extend(df['LMP'].tolist())
    else:
        print("file doesn't exist", file_name)
# print(len(all_G_values[0.0]))


# read load data on buses
def weights_AllDays():
    busfile = os.path.join(INDIR,'BusData.CSV')
    loadData = pd.read_csv(busfile)
    # calculate load weight
    loadData['Weight'] = loadData['Pd']/loadData['Pd'].sum()
    return loadData.set_index('Bus')['Weight'].to_dict()

# Function to load weights if have different weight for different day
def load_weights(d_value):
    weight_path = os.path.join(INDIR, 'BusData'+str(d_value)+'D.CSV')
    if os.path.exists(weight_path):
        weight_df = pd.read_csv(weight_path)
        # calculate load weight
        weight_df['Weight'] = weight_df['Pd']/weight_df['Pd'].sum()
        return weight_df.set_index('Bus')['Weight'].to_dict()
    else:
        raise FileNotFoundError(f"Weight file for d={d_value} not found.")


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

# Initialize a dictionary to store weighted average G values for each r
average_L_over_time = {}

# Iterate through each r value
for r in sorted(set(r for r, d in file_names_G.keys())):
    l_sum_over_t = None
    weight_sum_over_t = None
    count_d = 0
    # print("r: ", r)
    # Iterate through each (r, d) pair
    for (cur_r, cur_d), file_path in file_names_G.items():
        if cur_r == r:
            # print("r, d: ", cur_r, " ", cur_d," count d: ", count_d)
            if os.path.exists(file_path):
                df = pd.read_csv(file_path)  
                count_d += 1
                # print("count_d: ", count_d)
                # # Calculate average G values over periods t for the current (r, d) pair
                # avg_G_over_t = df.groupby('time')['LMP'].mean()
            
                # # Accumulate average G values over d
                # if avg_G_over_d is None:
                #    avg_G_over_d = avg_G_over_t
                # else:
                #    avg_G_over_d += avg_G_over_t
                
                # # Load weights for the current d value
                # # weights = load_weights(d)
                # weights = weights_AllDays()
                # df['Weight'] = df['busID'].map(weights)
                
                df['Weight'] = 1
                
                # Calculate weighted sum of L values over periods t for the current (r, d) pair
                weighted_L = df['LMP'] * df['Weight']
                sum_weighted_L = df.groupby('time')['Weight'].sum()
                weighted_avg_L_over_t = df.groupby('time').apply(lambda x: (x['LMP'] * x['Weight']).sum() / x['Weight'].sum())

                # print(weighted_avg_L_over_t)
                
                if l_sum_over_t is None:
                    # print("first day")
                    l_sum_over_t = weighted_avg_L_over_t
                    weight_sum_over_t = sum_weighted_L
                else:
                    l_sum_over_t += weighted_avg_L_over_t
                    weight_sum_over_t += sum_weighted_L
    
    # Average over d
    if l_sum_over_t is not None:
        avg_l_over_t = l_sum_over_t / count_d
        # Store the weighted averages in the dictionary
        average_L_over_time[r] = avg_l_over_t

    
# Plot histograms for each r value
plt.figure(figsize=(6, 10))
num_rows = int(np.ceil(len(all_G_values)))
for i, r in enumerate(all_G_values):
    plt.subplot(num_rows, 1, i + 1)
    mean_value = np.mean(average_L_over_time[r])
    # print(mean_value)
    plt.hist(all_G_values[r], bins=40, alpha=0.7, label=f'R = {r}', color='blue')
    plt.axvline(mean_value, color='r',linestyle='-',linewidth=3)
    plt.text(mean_value, -0.05, f'{mean_value:.2f}', va='top', ha='center', color='r')
    plt.xlabel('LMP ($/MWh)', fontname='Arial')
    plt.ylabel('Frequency', fontname='Arial')
    plt.legend()
    plt.grid(True)
    plt.xlim(0, 22)
    plt.ylim(0,230000)

plt.figtext(0.5, 0.0, f'Days {Day}', ha='center', va='bottom', fontsize=10, color='gray')
plt.tight_layout()
# plt.show()

# Plotting weighted average G vs periods t for each r
plt.figure(figsize=(12, 8))

i = 0
for r, avg_data in average_L_over_time.items():
    plt.plot(avg_data.index, avg_data.values, linestyle=linestyle_tuple[i], linewidth=5, marker='o', label=f'$R_b$ = {r}')
    i += 1
plt.xlabel('time (hour)', fontsize=35, fontname='Arial')
plt.ylabel('average LMPs ($/MWh)', fontsize=35, fontname='Arial')
# plt.title('Average LMPs at different bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')
# plt.title('Load-weighted average LMPs at different bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')
plt.xticks(np.arange(0,25,4),fontsize=35, fontname='Arial') 
plt.yticks(fontsize=35, fontname='Arial')
plt.ylim(3,16)
# plt.ylim(13,21)
plt.legend(fontsize=25)
plt.grid(True)
# plt.figtext(0.5, 0.01, f'Days {Day}', ha='center', va='bottom', fontsize=20, color='gray')
plt.tight_layout()
plt.show()





