import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import openpyxl
import os
import sys

INDIR = sys.argv[1]

file = os.path.join(INDIR,'GenNor_M10.xlsx')
# file = os.path.join(INDIR,'GenNor_MS15.xlsx')
df = pd.read_excel(file, sheet_name='AVG', index_col=0)
# print(df.head())

BRatio = [0.0,0.002,0.01,0.025,0.05,0.075,0.1,0.13]

# Plotting

linestyle_tuple = [
    (0, ()),
    (0, (3, 5, 1, 5)),
    (0, (1, 1)),
    (5, (10, 3)),
    (0, (5, 1)),
    
    (0,(2,2)),
    (0, (5, 10)),
    (0, (3, 10, 1, 10)),
    (0, (3, 1, 1, 1))]

c_idx = 0
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12,8))
ax1.plot(df.index, df["[0,0]"], label="[0,0]", color="blue", linestyle=linestyle_tuple[c_idx], linewidth=5)
ax1.set_ylim(0.82,0.85)  # Set the limits to emphasize around 0.83
ax1.spines['bottom'].set_visible(False)  # Hide the bottom spine
ax1.legend(loc='upper left') 
plt.yticks(fontsize=20, fontname='Arial') 

c_idx += 1
ax2.plot(df.index, df["[0,1000]"], label="[0,1000]", color="orange", linestyle=linestyle_tuple[c_idx], linewidth=5)
ax2.set_ylim(0.45,0.48)  # Set the limits to emphasize around 0.45
ax2.spines['top'].set_visible(False)  # Hide the top spine
ax2.legend(loc='upper left')  # Add legend for the second line

ax2.set_xlabel('Bidding ratio', fontsize=25, fontname='Arial')
ax2.set_ylabel('Normalized average generation rate', fontsize=25, fontname='Arial')
ax2.yaxis.set_label_coords(-0.08, 1)  # Positioning y-axis label at the center
ax1.tick_params(axis='both', which='major', labelsize=20)  
ax2.tick_params(axis='both', which='major', labelsize=20)  
fig.suptitle('Generation of "zero-cost" generators at different bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')


# plt.figure(figsize=(12, 8))
# for column in df.iloc[:,[0,2]].columns:
#     print(column)
#     plt.plot(df.index, df[column], label=column, linestyle=linestyle_tuple[c_idx], linewidth=5)
#     c_idx += 1
    
    
# plt.xlabel('Bidding ratio', fontsize=25, fontname='Arial')
# plt.ylabel('Normalized average generation rate', fontsize=25, fontname='Arial')


# # # Add a horizontal line 
# # plt.axhline(y=0.83, color='black', linestyle='--', label='Baseline (0.83)')
# # plt.axhline(y=0.71, color='black', linestyle=':', label='Baseline (0.71)')
# # plt.axhline(y=0.445, color='black', linestyle='-.', label='Baseline (0.445)')
    
# plt.title('Generation of "zero-cost" generators at different bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')

# plt.legend()
plt.xticks(BRatio)
# plt.ylim(0,1)
plt.xticks(fontsize=20, fontname='Arial') 
plt.yticks(fontsize=20, fontname='Arial') 
plt.xticks(rotation=60)  # Rotate x-axis labels for better readability
# Format x-axis tick labels to display float values with minimal precision
plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.3g}'))
plt.tight_layout()
plt.show()
