import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

INDIR = sys.argv[1]
BRatio = [0.0,0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.3]

set1 = 0
Setname = ['WZLT', 'WHLT', 'WFHLT', 'WTLT', 'WFLT', 'WTTL', 'WTTLT']
Setmeaning = [[0,0], [0,100], [0,500], [0,1000], [0,5000], [0,10000],[10000,10000]]
# Setname = ['OZLT', 'OHLT', 'OETL', 'OTLT', 'OWTL', 'OTTL','TTLT','ETTL']
# Setmeaning = [[0,0], [0,100], [0,800], [0,1000], [0,2000], [0,10000], [10000,10000], [8000,20000]]

# Create list of file names for two sets (LMPs and GenCost)
file_names_G = [None]*len(BRatio)
file_names_C = [None]*len(BRatio)
for i in range(len(BRatio)):
    BRnum = str(BRatio[i])
    file_names_G[i] = os.path.join(INDIR,'GLPAllbus_24h_' + Setname[set1] + BRnum + '.CSV')
    file_names_C[i] = os.path.join(INDIR,'Loss_24h_' + Setname[set1] + BRnum + '.xlsx')

# read load data on buses
busfile = os.path.join(INDIR,'BusData.CSV')
loadData = pd.read_csv(busfile)
# calculate load weight
loadData['Weight'] = loadData['Pd']/loadData['Pd'].sum()
loadWeight = loadData['Weight'].values

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

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# Read and plot data from set G in the first subplot
BRidx = 0
LMPBR = np.zeros(len(BRatio))
LMPPeak = np.zeros(len(BRatio))
LMPVar = np.zeros((len(BRatio),24))
LMPVarAVG = np.zeros(len(BRatio))
for file_G in file_names_G:
    data_G = pd.read_csv(file_G)
    buses = data_G.busID.unique()
    times = data_G.time.unique()
    LMPs = np.zeros((len(buses),24))
    for busi in range(len(buses)):
        databus = data_G.loc[data_G['busID']==buses[busi]]
        LMPs[busi,:] = databus.LMP
    LMPA = np.average(LMPs, axis=0, weights=loadWeight)
    LMPBR[BRidx] = np.average(LMPA)
    LMPPeak[BRidx] = np.max(LMPA)
    for ts in range(24):
        LMPVar[BRidx,ts] = np.sum((LMPs[:,ts]-LMPA[ts])**2,axis=0)
    LMPVarAVG = np.average(LMPVar,axis=1)
    axs[0].plot(times, LMPA, label=str(BRatio[BRidx]), linestyle=linestyle_tuple[BRidx], linewidth=5)
    BRidx += 1

# Set labels and title for the first subplot
axs[0].set_xlabel('time (hour)', fontsize=20, fontname='Arial')
axs[0].set_ylabel('LMPs ($/MWh)', fontsize=20, fontname='Arial')
axs[0].tick_params(axis='x', labelsize=15)
axs[0].tick_params(axis='y', labelsize=15)
# axs[0].set_xlim(1,24)

BRidx = 0
AVGCostBR = np.zeros(len(BRatio))
# Read and plot data from set C in the second subplot
for file_C in file_names_C:
    data_C = pd.read_excel(file_C, sheet_name='Loss', nrows=24)
    axs[1].plot(data_C['time'], data_C['GenCost']/data_C['GenTot'], label=str(BRatio[BRidx]) ,linestyle=linestyle_tuple[BRidx], linewidth=5)
    AVGCostBR[BRidx] = np.average(data_C['GenCost']/data_C['GenTot'])
    BRidx += 1

# Set labels and title for the second subplot
axs[1].set_xlabel('time (hour)', fontsize=20, fontname='Arial')
axs[1].set_ylabel('Average generation cost ($/MWh)', fontsize=20, fontname='Arial')
axs[1].tick_params(axis='x', labelsize=15)
axs[1].tick_params(axis='y', labelsize=15)
axs[1].legend(fontsize=12)
# axs[1].legend(bbox_to_anchor=(0.5, 1.15), borderaxespad=0, fancybox=True, shadow=True, ncol=4, fontsize=12)
# axs[1].set_xlim(1,24)
CostData = pd.DataFrame({"Bidding Ratio": BRatio, "Peak LMPs": LMPPeak, "AVG LMPs": LMPBR, "Varibility LMPs": LMPVarAVG, "AVG GenCost": AVGCostBR})
CostDataFile = os.path.join(INDIR,'CostData_24h_' + Setname[set1] + '.CSV')
CostData.to_csv(CostDataFile, index=False)
print("AVG LMPs: ", LMPBR)
print("Peak LMPs: ", LMPPeak)
print("AVG Cost: ", AVGCostBR)
print("Variability in LMPs: ", LMPVarAVG)
# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()
