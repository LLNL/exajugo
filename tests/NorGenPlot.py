import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import openpyxl
import os
import sys

INDIR = sys.argv[1]

file1 = os.path.join(INDIR,'GenNor_M10.xlsx')
file2 = os.path.join(INDIR,'GenNorAll.xlsx')
df1 = pd.read_excel(file1, sheet_name='AVG', index_col=0)
df2 = pd.read_excel(file2, sheet_name='AVG', index_col=0)
# print(df1.head())

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

# Generate a list of unique colors
colors = plt.cm.tab20(np.linspace(0, 1, 15))

c_idx = 0

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12,8))
ax1.plot(df1.index, df1["[0,0]"], label="Zero", color=colors[c_idx], linestyle=linestyle_tuple[c_idx], linewidth=5)
ax1.set_ylim(0.83,0.85)  # Set the limits to emphasize around 0.83
ax1.spines['bottom'].set_visible(False)  # Hide the bottom spine
# ax1.legend(loc='upper left') 
plt.yticks(fontsize=20, fontname='Arial') 


# fig, ax1 = plt.subplots(figsize=(12,8))

# line1, = ax1.plot(df1.index, df1["[0,0]"], label="Zero", color=colors[c_idx], linestyle=linestyle_tuple[c_idx], linewidth=5)
# ax1.set_ylim(0.83,0.85)  # Set the limits to emphasize around 0.83
# ax1.set_xlabel('Bidding ratio', fontsize=20, fontname='Arial')
# ax1.set_ylabel('Normalized average generation rate, zero-cost', fontsize=20, fontname='Arial',color=colors[c_idx])
# ax1.tick_params(axis='y', labelcolor=colors[c_idx],labelsize=20)
# # ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.3g}'))  # Format x-axis tick labels
# ax1.set_xticks(BRatio)
# ax1.tick_params(axis='x', rotation=60,labelsize=20)  # Rotate x-axis labels by 60 degrees
# ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.3g}'))  # Format x-axis tick labels

c_idx += 5
ax2.plot(df2.index, df2["[0,0]"], label="Nonzero", color=colors[c_idx], linestyle=linestyle_tuple[c_idx], linewidth=5)
ax2.set_ylim(0.17,0.18)  # Set the limits to emphasize around 
ax2.spines['top'].set_visible(False)  # Hide the top spine
# ax2.legend(loc='upper left')  # Add legend for the second line

ax2.set_xlabel('Bidding ratio', fontsize=25, fontname='Arial')
ax2.set_ylabel('Normalized average generation rate', fontsize=25, fontname='Arial')
ax2.yaxis.set_label_coords(-0.08, 1)  # Positioning y-axis label at the center
ax1.tick_params(axis='both', which='major', labelsize=20)  
ax2.tick_params(axis='both', which='major', labelsize=20)  
# ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.3g}'))
fig.suptitle('Generation of generators at different bidding ratio, 2000-bus setting', fontsize=25, fontname='Arial')


# ax2 = ax1.twinx()  # Creates a new y-axis that shares the same x-axis
# line2, = ax2.plot(df2.index, df2["[0,0]"], label="Nonzero", color=colors[c_idx], linestyle=linestyle_tuple[c_idx], linewidth=5)
# ax2.set_ylim(0.17,0.18)  # Set the limits to emphasize around
# ax2.set_ylabel('Normalized average generation rate, nonzero-cost', fontsize=20, fontname='Arial',color=colors[c_idx])
# ax2.tick_params(axis='y', labelcolor=colors[c_idx],labelsize=20)

# # lines = [line1, line2]
# # labels = [line.get_label() for line in lines]
# # ax1.legend(lines, labels, loc='upper left')


# plt.legend()
# plt.xticks(BRatio)
# # plt.ylim(0,1)
# plt.xticks(fontsize=20, fontname='Arial') 
# plt.yticks(fontsize=20, fontname='Arial') 
# plt.xticks(rotation=60)  # Rotate x-axis labels for better readability
# # Format x-axis tick labels to display float values with minimal precision
# plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.3g}'))
plt.tight_layout()
plt.show()
