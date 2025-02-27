import sys
import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

INDIR = sys.argv[1]
d = 1
Day = [1,35,80,93,355]
BRatio = np.array([0.0,0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3])



LMPA = np.zeros((len(BRatio),24))

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

# read load data on buses
busfile = os.path.join(INDIR,'BusData'+str(Day[d])+'D.CSV')
loadData = pd.read_csv(busfile)
# calculate load weight
loadData['Weight'] = loadData['Pd']/loadData['Pd'].sum()
loadWeight = loadData['Weight'].values

# read plant results, calculate AVG (load-weighted) LMPs, plot
for i in range(len(BRatio)):
    BRnum = str(BRatio[i])
    filename = 'GLPAllbus_24h_'+'TNOZ'+str(Day[d])+'D'+str(BRatio[i])+'.CSV'
    file = os.path.join(INDIR,filename)
    if os.path.exists(file):
        data = pd.read_csv(file)
        buses = data.busID.unique()

        LMPs = np.zeros((len(buses),24))
        for busi in range(len(buses)):
            databus = data.loc[data['busID']==buses[busi]]
            LMPs[busi,:] = databus.LMP
        
        LMPA[i,:] = np.average(LMPs, axis=0, weights=loadWeight)
        # LMPA[i,:] = np.mean(LMPs, axis=0)
        plt.plot(range(24),LMPA[i,:],label=str(BRatio[i]), linestyle=linestyle_tuple[i],linewidth=5)


plt.ylabel("LMPs ($/MWh)",fontsize=25, fontname="Arial")
plt.xlabel("Time (hr)",fontsize=25, fontname="Arial")
plt.title(f'Average LMPs at different bidding ratios [Day {Day[d]}]', fontsize=25, fontname="Arial")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
# plt.xlim(1,24)
plt.legend()

plt.show()



        





