#!/usr/bin/env python
# coding: utf-8
import os, sys, glob
import numpy as np
import pandas as pd
import fnmatch


def string_cleaner(tt):
    #Clean strings
    #Remove unicode characters
    tt = ''.join(i for i in tt if ord(i)<128)

    #Replace commas with semicolons
    tt = tt.replace('","',';')
    tt = tt.replace(',',';')
    return tt

on_nersc = True
# Set to directory with all the VI files to merge
if on_nersc:
  merged_dir = os.environ['HOME']+'/SV/VI_files/SV0/QSO_Andes/output/'
else:
  merged_dir = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO_Andes/output/'  

tiles = ['68002'] 
nights = ['20200315']  

# Read in list of files in merged directory
all_files = os.listdir(merged_dir)
merged_files=[]
pattern = "desi*"+tiles[0]+"*merged.csv"
for entry in all_files:
  if fnmatch.fnmatch(entry, pattern):
    merged_files.append(entry)
print(merged_files)


combined_file = merged_dir+"desi-vi_QSO_reinspection_"+tiles[0]+"_merged_all.csv"
print(combined_file)

# Read the first file in to vi to set up vi
print(merged_files[0])
vimerged = pd.read_csv(merged_dir + merged_files[0], delimiter = ",", engine='python')
# Read in the rest of the files and append them to vi
for i in range(1,len(merged_files)):
    print(merged_files[i])
    vi2 = pd.read_csv(merged_dir + merged_files[i], delimiter = ",", engine='python')
    vimerged = vimerged.append(vi2, ignore_index=True)


# Get rid of evil characters
vimerged['all VI comments'] = vimerged['all VI comments'].apply(string_cleaner)
vimerged['merger comment'] = vimerged['merger comment'].apply(string_cleaner)

#for i in np.arange(len(vimerged['TARGETID'])):
#  print(vimerged.loc[i]['all VI comments'])
#print(vimerged[['Redrock z', 'best z', 'best class']])

# Print to a combined file
vimerged[['TARGETID','Redrock_z', 'best z', 'best quality', 'Redrock_spectype', 'best spectype', 'all VI issues', 'all VI comments', 'merger comment','N_VI','DELTACHI2', 'ZWARN', 'ZERR','FIBER','FLUX_G', 'FLUX_R', 'FLUX_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MW_TRANSMISSION_G','MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z']].to_csv(combined_file,index=False)

vitest = pd.read_csv(combined_file)
print(vitest)
