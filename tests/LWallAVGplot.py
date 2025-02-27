import sys
import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

INDIR = sys.argv[1]
BRatio = [0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.3]
# BRatio = [0.001,0.01,0.025,0.05,0.1]

set1 = 0
Setname = ['WZLT', 'WHLT', 'WFHLT', 'WTLT', 'WFLT', 'WTTL', 'WTTLT']
Setmeaning = [[0,0], [0,100], [0,500], [0,1000], [0,5000], [0,10000],[10000,10000]]
# Setname = ['OZLT', 'OHLT', 'OETL', 'OTLT', 'OWTL', 'OTTL','TTLT','ETTL']
# Setmeaning = [[0,0], [0,100], [0,800], [0,1000], [0,2000], [0,10000], [10000,10000], [8000,20000]]

LMPA = np.zeros((len(BRatio),24))

linestyle_tuple = [
    (0, (1, 1)),
    (0, (3, 5, 1, 5)),
    (5, (10, 3)),
    (0, (5, 1)),
    (0, ()),
    
    (0, (5, 10)),
    (0, (3, 10, 1, 10)),
    (0, (3, 1, 1, 1))]

# read load data on buses
busfile = os.path.join(INDIR,'BusData.CSV')
loadData = pd.read_csv(busfile)
# calculate load weight
loadData['Weight'] = loadData['Pd']/loadData['Pd'].sum()
loadWeight = loadData['Weight'].values

# read plant results, calculate AVG (load-weighted) LMPs, plot
for i in range(len(BRatio)):
    BRnum = str(BRatio[i])
    filename = 'GLPAllbus_24h_' + Setname[set1] + BRnum + '.CSV'
    file = os.path.join(INDIR,filename)
    data = pd.read_csv(file)
    buses = data.busID.unique()

    LMPs = np.zeros((len(buses),24))
    for busi in range(len(buses)):
        databus = data.loc[data['busID']==buses[busi]]
        LMPs[busi,:] = databus.LMP
        
    LMPA[i,:] = np.average(LMPs, axis=0, weights=loadWeight)
    # LMPA[i,:] = np.mean(LMPs, axis=0)
    plt.plot(range(24),LMPA[i,:],label=str(BRatio[i]), linestyle=linestyle_tuple[i],linewidth=5)


plt.ylabel("LMPs ($/MWh)",fontsize=30, fontname="Arial",)
plt.xlabel("Time (hr)",fontsize=30, fontname="Arial",)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
# plt.xlim(1,24)
plt.legend(bbox_to_anchor=(0.5, 1.15),loc='upper center', borderaxespad=0,fancybox=True, shadow=True, ncol=4, fontsize=20)

plt.show()



        





