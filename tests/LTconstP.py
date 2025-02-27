import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

INDIR = sys.argv[1]


# BRatio = [0.0,0.002,0.01,0.025,0.05,0.075,0.1,0.13]
BRatio = [0.0,0.13]
Day = 224
set1 = 0
set2 = 2
phase = 1


Setname = ['M10OZ', 'M10HOZ', 'M10TOZ', 'M10TTOZ']
Setmeaning = [[0,0], [0,100], [0,1000], [0,10000]]

# Create list of file names for two sets 
file_names_A = [None]*len(BRatio)
file_names_B = [None]*len(BRatio)
for i in range(len(BRatio)):
    # print("file i ",i)
    BRnum = str(BRatio[i])
    file_names_A[i] = os.path.join(INDIR,'Loss_24h_' + Setname[set1] + str(Day) + 'D' + BRnum + '.xlsx')
    file_names_B[i] = os.path.join(INDIR,'Loss_24h_' + Setname[set2] + str(Day) + 'D' + BRnum + '.xlsx')

# BRatio = [0.0,0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.3]
# set1 = 0
# set2 = 4
# phase = 1


# Setname = ['WZLT', 'WHLT', 'WFHLT', 'WTLT', 'WFLT', 'WTTL', 'WTTLT']
# Setmeaning = [[0,0], [0,100], [0,500], [0,1000], [0,5000], [0,10000],[10000,10000]]
# # Setname = ['OZLT', 'OHLT', 'OETL', 'OTLT', 'OWTL', 'OTTL','TTLT','ETTL']
# # Setmeaning = [[0,0], [0,100], [0,800], [0,1000], [0,2000], [0,10000], [10000,10000], [8000,20000]]

# # Create list of file names for two sets 
# file_names_A = [None]*len(BRatio)
# file_names_B = [None]*len(BRatio)
# for i in range(len(BRatio)):
#     BRnum = str(BRatio[i])
#     file_names_A[i] = os.path.join(INDIR,'Loss_24h_' + Setname[set1] + BRnum + '.xlsx')
#     file_names_B[i] = os.path.join(INDIR,'Loss_24h_' + Setname[set2] + BRnum + '.xlsx')

phase1 = 1
phase2 = 2
# Create a separate figure for each sheet (L and T)
for sheet_name in ['Lactive', 'Tactive']:
    # Create a new figure for each sheet
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    # # Initialize an empty list to store all x-values
    # all_x_values = []
    
    for i, (file_name_A, file_name_B) in enumerate(zip(file_names_A, file_names_B)):
        # print("i ", i)
        # Read specific sheets from the Excel files in set A and set B
        sheet_A = pd.read_excel(file_name_A, sheet_name=sheet_name)
        sheet_B = pd.read_excel(file_name_B, sheet_name=sheet_name)

        # Filter data by the specified column
        # combine phases 1 and 2
        # data_A2 = sheet_A
        # data_B2 = sheet_B
        data_A2 = sheet_A[sheet_A['phase']==phase]
        data_B2 = sheet_B[sheet_B['phase']==phase]


        # Extract x and y coordinates from the filtered data
        if sheet_name == 'Lactive':
            # # Append unique x-values to the list
            # all_x_values.extend(data_A2['line'].unique())
            # all_x_values.extend(data_B2['line'].unique())
            x_A2 = data_A2['line'].values
            x_B2 = data_B2['line'].values
        else:
            # all_x_values.extend(data_A2['trans'].unique())
            # all_x_values.extend(data_B2['trans'].unique())
            x_A2 = data_A2['trans'].values
            x_B2 = data_B2['trans'].values
        y_A2 = data_A2['status'].values
        y_B2 = data_B2['status'].values
            
        # # Get unique x-values from the list
        # unique_x_values = sorted(set(all_x_values))
        
        # Calculate subplot indices
        row = i // 2
        col = i % 2
        # print(row," ",col)
        # Plot the data from set A and set B in the current subplot
        axs[i].scatter(x_A2, y_A2, alpha=0.4, marker='o', label=Setmeaning[set1], color='red')
        axs[i].scatter(x_B2, y_B2, alpha=0.8, marker='x', label=Setmeaning[set2], color='blue')
        # axs[row, col].scatter(x_A2, y_A2, alpha=0.8, marker='x', label=Setmeaning[set1], color='blue')
        # axs[i].set_xticks(unique_x_values)
        # # Rotate the x-labels to prevent overlapping
        # axs[i].set_xticklabels(unique_x_values, rotation=45)
        
        # axs[i].set_xticks(unique_x_values)
        # # Rotate the x-labels to prevent overlapping
        # axs[i].set_xticklabels(unique_x_values, rotation=45)
        
        # Set labels and title for the subplot
        axs[i].set_xlabel('line #' if sheet_name == 'Lactive' else 'transformer #', fontname='Arial',fontsize=20)
        axs[i].set_ylabel('time (hour)', fontname='Arial',fontsize=20)
        axs[i].set_title(f'R= {BRatio[i]}',fontname='Arial',fontsize=20)
        axs[i].set_ylim(0,25)
        axs[i].set_yticks([0,4,8,12,16,20,24])
        axs[i].tick_params(axis='both', which='major', labelsize=15)
        # axs[row, col].legend()

    # Get handles and labels from the last subplot
    handles, labels = axs[0].get_legend_handles_labels()

    # Add a single legend outside of the subplots
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=2)
    plt.figtext(0.5, 0.01, f'Day {Day}', ha='center', va='bottom', fontsize=15, color='gray')
    # Set the title for the figure
    fig.suptitle(f'Filtered by {"phase"} = {phase}')

    # Adjust layout to prevent overlap and increase space between subplots
    plt.subplots_adjust(hspace=0.5)
    
    # # Adjust layout to prevent overlap
    # plt.tight_layout()

    # Show the plot for the current sheet
    plt.show()
