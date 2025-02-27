import sys
import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

INDIR = sys.argv[1]
BRatio = [0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.3]
# BRatio = [0.01]
# BRatio = [0.001,0.01,0.025,0.05,0.1]

set1 = 4
Setname = ['WZLT', 'WHLT', 'WFHLT', 'WTLT', 'WFLT', 'WTTL', 'WTTLT']
Setmeaning = [[0,0], [0,100], [0,500], [0,1000], [0,5000], [0,10000],[10000,10000]]
# Setname = ['OZLT', 'OHLT', 'OETL', 'OTLT', 'OWTL', 'OTTL','TTLT','ETTL']
# Setmeaning = [[0,0], [0,100], [0,800], [0,1000], [0,2000], [0,10000], [10000,10000], [8000,20000]]

LMPA = np.zeros((len(BRatio),24))
PPA = np.zeros((len(BRatio),24))

linestyle_tuple = [
    (0, (3, 5, 1, 5)),
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

# read plant results, calculate AVG (load-weighted) LMPs, plot
for i in range(len(BRatio)):
    BRnum = str(BRatio[i])
    filename = 'Plantbus_24h_' + Setname[set1] + BRnum + '.CSV'
    file = os.path.join(INDIR,filename)
    data = pd.read_csv(file)
    buses = data.busID.unique()
    # print(buses)
    # calculate load weight
    loadbus = loadData.loc[loadData['Bus'].isin(buses)]
    loadbus.loc[:,'Weight'] = loadbus['Pd']/loadbus['Pd'].sum()
    loadWeight = loadbus['Weight'].values
    # print(loadbus.head(30))
    # print(loadbus.shape)
    LMPs = np.zeros((len(buses),24))
    pplant = np.zeros((len(buses),24))
    for busi in range(len(buses)):
        databus = data.loc[data['busID']==buses[busi]]
        LMPs[busi,:] = databus.LMP
        pplant[busi,:] = databus.plantpower
    LMPA[i,:] = np.average(LMPs, axis=0, weights=loadWeight)
    PPA[i,:] = np.average(pplant, axis=0, weights=loadWeight)
    # plt.plot(range(24),PPA[i,:],label='R=%.4f' %BRatio[i], linestyle=linestyle_tuple[i],linewidth=5)
    plt.plot(range(24),LMPA[i,:],label='R=%.4f' %BRatio[i], linestyle=linestyle_tuple[i],linewidth=5)


plt.ylabel("LMPs ($/MWh)",fontsize=30, fontname="Arial",)
# plt.ylabel("Plant power (MW)",fontsize=30, fontname="Arial",)
plt.xlabel("Time (hr)",fontsize=30, fontname="Arial",)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(bbox_to_anchor=(0.5, 1.15),loc='upper center', borderaxespad=0,fancybox=True, shadow=True, ncol=4, fontsize=20)

plt.show()



        





