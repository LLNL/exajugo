import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy import stats

INDIR = sys.argv[1]
# Day = [1,35]
BRatio = [0.0,0.001,0.01,0.025,0.05,0.1]
Day = [1,35,176,185,195,223,268,355]
# BRatio = [0.0,0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3]

# Calculate mean and std for each (d, r) pair
mean_values = []
sem_values = []
std_values = []
d_values = []

# Create list of file names and read data 
for d in Day:
    for r in BRatio:
        file = os.path.join(INDIR,'Loss_24h_'+'NOZ'+str(d)+'D'+str(r)+'.xlsx')

        if os.path.exists(file):
            df = pd.read_excel(file, sheet_name='Loss', nrows=24)  # Use first column as index (periods)
            df['GencostAvg'] = df['GenCost']/df['GenTot']
            t_values = df['time']
            G_values = df['GencostAvg'] 
            
            # Calculate mean and std for G
            mean = np.mean(G_values)
            sem = stats.sem(G_values)
            std = np.std(G_values)
            
            # Append to lists
            mean_values.append(mean)
            sem_values.append(sem)
            std_values.append(std)
            d_values.append(d)

# Convert to numpy arrays for easier manipulation
mean_values = np.array(mean_values)
sem_values = np.array(sem_values)
std_values = np.array(std_values)
d_values = np.array(d_values)
    
# Initialize arrays to store mean and overall mean for each d
mean_values_per_d = []
overall_mean_values = []

# Calculate mean for each d across all r values
for d_val in np.unique(d_values):
    indices_d = np.where(d_values == d_val)[0]
    mean_values_d = mean_values[indices_d]
    mean_values_per_d.append(mean_values_d)

# Convert lists to numpy arrays
mean_values_per_d = np.array(mean_values_per_d)
# print(mean_values_per_d)
overall_mean_values = np.mean(mean_values_per_d,axis=0)
overall_mean_values = np.array(overall_mean_values)   
# print("overall: ", overall_mean_values)
overall_sem = stats.sem(mean_values_per_d)             
        
linestyle_tuple = [
    # (0, (3, 5, 1, 5)),
    (0, (1, 1)),
    (0, (3, 5, 1, 5)),
    (5, (10, 3)),
    (0, (5, 1)),
    (0, ()),
    (0,(2,2)),
    (0, (5, 10)),
    (0, (3, 10, 1, 10)),
    (0, (3, 1, 1, 1))]

# plot mean lines, standard deviation bands, and range bands
plt.figure(figsize=(10, 6))
colors = plt.cm.tab20(np.linspace(0, 1, len(Day)))
# Plot lines for each d with shaded regions for standard deviation
i = 0
for d_value in np.unique(d_values):
    indices = np.where(d_values == d_value)[0]
    plt.plot(BRatio, mean_values[indices], label=f'd = {d_value}',color=colors[i],linestyle=linestyle_tuple[i],linewidth=3)
    plt.errorbar(BRatio, mean_values[indices], yerr=1.96 * sem_values[indices],color=colors[i])
    plt.fill_between(BRatio, mean_values[indices] - std_values[indices], mean_values[indices] + std_values[indices], alpha=0.2,color=colors[i])
    i += 1
# Plot overall mean across all d values

plt.plot(BRatio, overall_mean_values, linewidth=2.5, color='red', label='Overall Mean')
plt.fill_between(BRatio, overall_mean_values - 1.96 * overall_sem, 
                 overall_mean_values + 1.96 * overall_sem, color='red', alpha=0.2)

# Customize plot
plt.xlabel('Bidding ratio', fontsize=25, fontname='Arial')
plt.ylabel('Average gen cost ($/MWh)', fontsize=25, fontname='Arial')
plt.title('Average generation cost at differnt bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')
plt.legend()
plt.xticks(BRatio)


plt.xticks(fontsize=20, fontname='Arial') 
plt.yticks(fontsize=20, fontname='Arial') 
plt.xticks(rotation=60)  # Rotate x-axis labels for better readability
# Format x-axis tick labels to display float values with minimal precision
plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.3g}'))
plt.tight_layout()
plt.show()



                
        
    


