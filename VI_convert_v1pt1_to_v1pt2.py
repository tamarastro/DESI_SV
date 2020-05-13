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

def comment_compressor(tt):
	# Get rid of superfluous --|--|... and all |
	tt = tt.replace("--|","")
	tt = tt.replace("|"," ")
	# Remove leading spaces
	tt = tt.strip()
	return tt

directory_truth_tables = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO/truth_tables/'
list_of_truth_tables_v1pt1 = ['truth_table_BGS_v1.1.csv','truth_table_ELG_v1.1.csv','truth_table_LRG_v1.1.csv']
list_of_truth_tables_v1pt2 = ['truth_table_BGS_v1.2.csv','truth_table_ELG_v1.2.csv','truth_table_LRG_v1.2.csv']

# Convert the version 1.1 files to version 1.2
for i,table in enumerate(list_of_truth_tables_v1pt1):
	vi = pd.read_csv(directory_truth_tables+list_of_truth_tables_v1pt1[i], delimiter = ",", engine='python')

	# Change the column names to v1.2
	vi = vi.rename(columns={"best flag": "best quality"})
	vi = vi.rename(columns={"all VI issue": "all VI issues"})

	# Get rid of superfluous --|-- and |
	vi.loc[vi['all VI comments']!='--', 'all VI comments'] = vi.loc[vi['all VI comments']!='--', 'all VI comments'].apply(comment_compressor)
	vi.loc[vi['all VI issues']!='--', 'all VI issues'] = vi.loc[vi['all VI issues']!='--', 'all VI issues'].apply(comment_compressor)
	
	# Get rid of duplicate issue flags	
	vi.loc[vi['all VI issues']!='--', 'all VI issues'] = vi.loc[vi['all VI issues']!='--', 'all VI issues'].transform(lambda x: ''.join(set(list(x))-set('-')))
	vi.loc[vi['all VI issues']!='--', 'all VI issues'] = vi.loc[vi['all VI issues']!='--', 'all VI issues'].transform(lambda x: x.replace(" ",""))
	
	# Get rid of unicode and commas
	vi['all VI comments']=vi['all VI comments'].apply(string_cleaner)
	vi['merger comment']=vi['merger comment'].apply(string_cleaner)

	# Save to file without the pandas index at the front and with slightly rearranged columns
	vi[['TARGETID','Redrock z', 'best z', 'best quality', 'Redrock spectype', 'best spectype', 'all VI issues', 
	'all VI comments', 'merger comment','N_VI','DELTACHI2', 'ZWARN', 'ZERR',
	'FIBER','FLUX_G', 'FLUX_R', 'FLUX_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 
	'MW_TRANSMISSION_G','MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z'
	]].to_csv(directory_truth_tables+list_of_truth_tables_v1pt2[i],index=False)
