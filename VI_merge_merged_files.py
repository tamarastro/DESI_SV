#!/usr/bin/env python
# coding: utf-8
import os, sys, glob
import numpy as np
import pandas as pd
import fnmatch

merged_dir = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO/output/' #dir with vi merged output
tiles = ['68002']
nights = ['20200315']
#subset = "_3_" 

# Read in list of files in merged directory
all_files = os.listdir(merged_dir)
merged_files=[]
pattern = "desi*merged.csv"
for entry in all_files:
  if fnmatch.fnmatch(entry, pattern):
    merged_files.append(entry)
print(merged_files)


combined_file = merged_dir+"desi-vi_SV0_QSO_tile"+tiles[0]+"_night"+nights[0]+"_merged_all.csv"
print(combined_file)

# Read the first file in to vi to set up vi
print(merged_files[0])
vimerged = pd.read_csv(merged_dir + merged_files[0], delimiter = ",", engine='python')
# Read in the rest of the files and append them to vi
for i in range(1,len(merged_files)):
    print(merged_files[i])
    vi2 = pd.read_csv(merged_dir + merged_files[i], delimiter = ",", engine='python')
    vimerged = vimerged.append(vi2, ignore_index=True)

#for i in np.arange(len(vimerged['TargetID'])):
#  print(vimerged.loc[i]['all VI comments'])
#print(vimerged[['Redrock z', 'best z', 'best class']])

# Print to a combined file
vimerged[['TargetID','Redrock z', 'best z', 'best class', 'Redrock spectype', 'best spectype', 'best issue', 'all VI comments', 'merger comment', 'N_VI', 'DELTACHI2','FIBER','FLUX_G','FLUX_R','FLUX_Z','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z']].to_csv(combined_file,index=False)

vitest = pd.read_csv(combined_file)
print(vitest)
