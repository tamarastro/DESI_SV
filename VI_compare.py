#!/usr/bin/env python
# coding: utf-8
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import fnmatch
import argparse

def plot_scatter_colors(x,y,col,clim,xlabel,ylabel,collabel,plot_name,xlim=0,ylim=0,psize=8.0,iextra=[9999],annotation=''):
  plt.figure(dpi=200)
  # Plot the VI disagreements using crosses
  if iextra[0] != 9999:
      plt.scatter(x[iextra],y[iextra], marker='x',s=18.0,
              c=col[iextra], cmap=plt.cm.get_cmap('jet_r'))
  plt.scatter(x,y, s=psize,
              c=col, cmap=plt.cm.get_cmap('jet_r'),alpha=0.9)
  plt.colorbar(ticks=range(5), label=collabel)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if xlim !=0: plt.xlim(xlim)
  if ylim !=0: plt.ylim(ylim)
  plt.annotate(annotation,xy=(0.05,0.9),xycoords='axes fraction')
  plt.clim(clim[0],clim[1])
  plt.savefig(plot_name,bbox_inches='tight')
  plt.close()

def plot_histogram(x,bin_edges,xlabel,ylabel,plot_name,comp=[999],annotation='',title=''):
  plt.figure(dpi=200)
  plt.hist(x, bins = bin_edges,rwidth=0.85)
  if comp[0] != 999:
    plt.hist(comp, bins = bin_edges,rwidth=0.25)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.annotate(annotation,xy=(0.05,0.9),xycoords='axes fraction')
  plt.savefig(plot_name,bbox_inches='tight')
  plt.close()

def compare_vi(VI_dir,VI_file,VI_truth_file,output_dir):
  VI_base = VI_file[:-4]  # Get rid of .csv at the end of the filename so we can add .png etc... for figures and output. 
  # Compare different visual inspections
  t = pd.read_csv(VI_truth_file) #truth (merged) file or file you want to compare to
  v = pd.read_csv(VI_dir+VI_file,delimiter = " , ", engine='python') #NOTE, prospect outputs data with spaces between commas, so you need the weird delimeter for prospect output, but not for normal csv output like we have done for the merged files
  # Change the column name to TARGETID to match standards elsewhere in DESI.
  v = v.rename(columns={"TargetID": "TARGETID"})

  # If the VI agreed with Redrock, insert the Redrock results into the VI column.
  v.loc[v['VI z']=='--', 'VI z'] = v.loc[v['VI z']=='--', 'Redrock z']
  v['VI z'] = v['VI z'].astype(float)
  v.loc[v['VI spectype']=='--', 'VI spectype'] = v.loc[v['VI spectype']=='--', 'Redrock spectype']

  # Gather the results by TARGETID, put them all in g, including only those items the VI'er looked at.
  g = v.merge(t,how='left',on='TARGETID')

  # Find the disagreements
  i_disagree_z = np.abs(g['VI z']-g['best z']) >=0.0033 
  #print( g[['TARGETID','best z','VI z']][i_disagree_z])
  #print(sum(i_disagree_z))

  i_disagree_class = abs(g['VI class']-g['best quality']) >= 2  
  #print( g[['TARGETID','best quality','VI class']][i_disagree_class])
  #print(sum(i_disagree_class))

  i_disagree_spectype = g['VI spectype'] != g['best spectype']
  #print( g[['TARGETID','best spectype','VI spectype']][i_disagree_spectype])
  #print(sum(i_disagree_spectype))

  i_good_bestz = g['best quality'] >=2.5  # Choose only those redshifts that were given a quality >=2.5 in the merged results
  i_good_VIz   = g['VI class']   >=3.0  # Choose only those redshifts that were given a quality >=3 by the VI
  i_good_z = i_good_bestz | i_good_VIz

  i_disagree_all = i_disagree_class | i_disagree_z | i_disagree_spectype
  i_disagree_good = i_disagree_z & i_good_z

  #----------------
  # PRINT TO SCREEN
  print('\nOut of %s objects, you disagreed with the merged results on %s redshifts, %s qualities, and %s spectypes.'%(len(g['TARGETID']),sum(i_disagree_z),sum(i_disagree_class),sum(i_disagree_spectype)))
  print('\nNumber of redshift disagreements on high-quality redshifts = %s.'%sum(i_disagree_good))
  print('\nTotal number of disagreements = %s.'%sum(i_disagree_all))
  print( g[['TARGETID','best z','VI z','best quality','VI class','best spectype','VI spectype']][i_disagree_all])

  conflicts = ''
  for target in g['TARGETID'][i_disagree_all]:
    conflicts = conflicts+', '+str(target)
  print('\nLook at your disagreements by opening the Prospect notebook viewer https://github.com/desihub/prospect/blob/master/doc/nb/Prospect_TARGETID.ipynb and cutting and pasting the following target list:')
  print(conflicts[2:])

  #-----------------------
  # WRITE TO FILE AND PLOT
  summaryfile = open(output_dir+VI_base+'_results.txt', "w") ####### PRINT OUTPUT TO FILE ########
  summaryfile.write('\nOut of %s objects, you disagreed with the merged results on %s redshifts, %s qualities, and %s spectypes.'%(len(g['TARGETID']),sum(i_disagree_z),sum(i_disagree_class),sum(i_disagree_spectype)))
  summaryfile.write('\nNumber of redshift disagreements on high-quality redshifts = %s.'%sum(i_disagree_good))
  summaryfile.write('\nTotal number of disagreements = %s.\n'%sum(i_disagree_all))
  #gd=g[i_disagree_all]
  g[['TARGETID','best z','VI z','best quality','VI class','best spectype','VI spectype']][i_disagree_all].to_string(summaryfile)
  summaryfile.write('\n\nLook at your disagreements by opening the Prospect notebook viewer https://github.com/desihub/prospect/blob/master/doc/nb/Prospect_TARGETID.ipynb and cutting and pasting the following target list:\n')
  summaryfile.write(conflicts[2:])
  summaryfile.close()
  # Plot VI z vs Merged z and colour by Quality
  plot_scatter_colors(g['best z'],g['VI z'],g['best quality'],[0.0,4.0],'Best z','VI z','Quality (4 is good)',output_dir+VI_base+'_z.png',iextra=i_disagree_good,annotation='# of z disagreements on high-quality z: %s'%sum(i_disagree_good))
  plot_histogram(g['best quality'],np.arange(11)/2.-0.25,'Quality','Number of cases',output_dir+VI_base+'_quality.png',comp=g['VI class'],annotation='Blue=Merged VI; Orange=Single VI')


#------------------------------------------------------------------
# Read in truth table (or any output file you'd like to compare to)
#VI_truth_file = '/global/cfs/cdirs/desi/sv/vi/TruthTables/truth_table_QSO_tile68002_night20200315.csv'
VI_truth_file  = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO/output/desi-vi_SV0_QSO_tile68002_night20200315_1_merged.csv'

# Define the VI file you want to compare, and the directory you find it in, and the output directory you want to send it to
VI_file    =  'desi-vi_SV0_QSO_tile68002_night20200315_1_TMD.csv'   # We remove .csv and add e.g. _z.png to print plots for this file
VI_dir     = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO/'
output_dir = VI_dir+'VI_check/'  # The directory in which we'll put the output. 
if not os.path.exists(output_dir):
  os.makedirs(output_dir)

# Set this to True if you want to do a batch job!  
# This will ignore the vi file you enter as VI_file above, but uses the directories and truth table.
Compare_Everything_In_Directory = False


# Note, the compare everything in directory option is fragile... requires you to set the pattern you want below.
if Compare_Everything_In_Directory:  
  # Find some VI files to check
  all_files = os.listdir(VI_dir)
  vi_files=[]
  # Choose a subset
  subset = '_1_'
  pattern = "desi*"+subset+"*.csv"
  for entry in all_files:
    if fnmatch.fnmatch(entry, pattern):
      vi_files.append(entry)
else:
  vi_files=[VI_file]


print(vi_files)

for file in vi_files:
  compare_vi(VI_dir,file,VI_truth_file,output_dir)
