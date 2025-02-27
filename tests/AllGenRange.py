import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import openpyxl
import os
import sys

DATADIR = sys.argv[1]

df = pd.read_excel(os.path.join(DATADIR,'GenAll.xlsx'))

ranges = ['Plb','Pub','Qlb','Qub','p0','q0']
# ranges = df.columns[3:9]
grouped = df.groupby('Bus',as_index=False).agg(
    {col: 'sum' for col in ranges} | {col: 'first' for col in df.columns if col not in ranges}
)


grouped.to_csv(os.path.join(DATADIR,'GenAllmz.csv'), index=False)
print("end of ")