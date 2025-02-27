import sys
import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

INDIR = sys.argv[1]
BRatio = [0.001,0.01,0.025,0.05,0.1,0.15,0.2,0.3]
# BRatio = [0.001,0.01,0.025,0.05,0.1,0.15,0.17]
# BRatio = [0.001,0.01,0.025,0.05,0.1]
# BRatio = [0.001, 0.01, 0.1]
LMPA = np.zeros((len(BRatio),24))
PPA = np.zeros((len(BRatio),24))

linestyle_tuple = [
    (0, (1, 1)),
    (0, (3, 5, 1, 5)),
    (5, (10, 3)),
    (0, (5, 1)),
    (0, ()),
    
    (0, (5, 10)),
    (0, (3, 10, 1, 10)),
    (0, (3, 1, 1, 1))]

for i in range(len(BRatio)):
    BRnum = str(BRatio[i])
    filename = 'Plantbus_24h_OZTTS' + BRnum + '.CSV'
    file = os.path.join(INDIR,filename)
    data = pd.read_csv(file)
    buses = data.busID.unique()
    LMPs = np.zeros((len(buses),24))
    pplant = np.zeros((len(buses),24))
    for busi in range(len(buses)):
        databus = data.loc[data['busID']==buses[busi]]
        LMPs[busi,:] = databus.LMP
        pplant[busi,:] = databus.plantpower
    LMPA[i,:] = np.mean(LMPs, axis=0)
    PPA[i,:] = np.mean(pplant, axis=0)
    plt.plot(range(24),PPA[i,:],label='R=%.4f' %BRatio[i], linestyle=linestyle_tuple[i],linewidth=5)
    # plt.plot(range(24),LMPA[i,:],label='R=%.4f' %BRatio[i], linestyle=linestyle_tuple[i],linewidth=5)

# i=2
# fig, ax1 = plt.subplots()
# color1 = 'tab:red'
# color2 = 'tab:blue'
# ax1.set_xlabel("Time (hr)",fontsize=16)
# ax1.set_ylabel("Plant power (MW)", fontsize=16, color=color1)
# ax1.plot(range(24),PPA[i,:], color=color1, linewidth=2.5)
# ax1.tick_params(axis='y', labelcolor=color1)
# ax1.set_title('Power and LMPs for R=%.4f' % BRatio[i])

# ax2 = ax1.twinx()
# ax2.set_ylabel("LMPs ($/MWh)",fontsize=16,color=color2)
# ax2.plot(range(24),LMPA[i,:], color=color2, linestyle=linestyle_tuple[1], linewidth=2.5)
# ax2.tick_params(axis='y', labelcolor=color2)
# fig.tight_layout

# LMPA = np.transpose(np.hstack((np.transpose(np.asmatrix(BRatio)),LMPA)))
# PPA = np.transpose(np.hstack((np.transpose(np.asmatrix(BRatio)),PPA)))
# fileLMP = os.path.join(INDIR,'LMPA.CSV')
# filePP = os.path.join(INDIR,'PPA.CSV')
# np.savetxt(fileLMP, LMPA, delimiter=",")
# np.savetxt(filePP, PPA, delimiter=",")

# plt.ylabel("LMPs ($/MWh)",fontsize=30, fontname="Arial",)
plt.ylabel("Plant power (MW)",fontsize=30, fontname="Arial",)
plt.xlabel("Time (hr)",fontsize=30, fontname="Arial",)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(bbox_to_anchor=(0.5, 1.15),loc='upper center', borderaxespad=0,fancybox=True, shadow=True, ncol=4, fontsize=20)

plt.show()



        





