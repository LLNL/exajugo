import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# old figure
INDIR = sys.argv[1]

Day = [1,35,80,93,355]
BRatio = [0.0,0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3]

# Dictionary to store data organized by BRatio values
dataGCost = {}

# Create list of file names and read data 
for r in BRatio:
    dfs = []
    r_key = (r,)
    for d in Day:
        file = os.path.join(INDIR,'Loss_24h_'+'NZ'+str(d)+'D'+str(r)+'.xlsx')
        
        if os.path.exists(file):
            df = pd.read_excel(file, sheet_name='Loss', nrows=24, index_col=0)  # Use first column as index (periods)
            df['GencostAvg'] = df['GenCost']/df['GenTot']
            dfs.extend(df['GencostAvg'].values)
    
    if dfs:
        dataGCost[r_key] = np.array(dfs)
                
        
# Calculate the mean, standard deviation, and range for each BRatio for each day
aggregated_data = {}
for r_key, g_values in dataGCost.items():
    r_value = r_key[0]  # Extract float value from tuple
    
    mean_g = np.mean(g_values)
    std_g = np.std(g_values)
    min_g = np.min(g_values)
    max_g = np.max(g_values)
    
    aggregated_data[r_key] = {
        'mean': mean_g,
        'std': std_g,
        'min': min_g,
        'max': max_g,
        'all_g_values': g_values 
    }

# plot mean lines, standard deviation bands, and range bands
plt.figure(figsize=(10, 6))

for r_key, stats in aggregated_data.items():
    r_value = r_key[0]  # Extract float value from tuple
    
    mean_g = stats['mean']
    std_g = stats['std']
    min_g = stats['min']
    max_g = stats['max']
    all_g_values = stats['all_g_values']
    
    # Plot data line with points 
    plt.plot([r_value], [mean_g], marker='o', markersize=8, color='blue', linewidth=5)
    # Error bars 
    plt.errorbar([r_value], [mean_g], yerr=std_g, fmt='', color='blue', capsize=5)
    # Fill between (blue shaded region with transparency)
    plt.fill_between([r_value], min_g, max_g, color='blue', alpha=0.3)
     # Scatter plot of all G values (gray)
    plt.scatter([r_value] * len(all_g_values), all_g_values, color='gray', alpha=0.5)  # Scatter plot of all G values

# Customize plot
plt.xlabel('Bidding ratio', fontsize=25, fontname='Arial')
plt.ylabel('Average gen cost ($/MWh)', fontsize=25, fontname='Arial')
plt.title('Average generation cost at differnt bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')

plt.xticks(BRatio)


plt.xticks(fontsize=20, fontname='Arial') 
plt.yticks(fontsize=20, fontname='Arial') 
plt.xticks(rotation=60)  # Rotate x-axis labels for better readability
# Format x-axis tick labels to display float values with minimal precision
plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.3g}'))
plt.tight_layout()
plt.show()



                
        
    


