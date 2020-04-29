#!/usr/bin/env python
# coding: utf-8
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

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

# Read in truth table (or any output file you'd like to compare to)
#VI_truth = '/global/cfs/cdirs/desi/sv/vi/TruthTables/truth_table_QSO_tile68002_night20200315.csv'
VI_dir     = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO/'
VI_merged  = 'output/desi-vi_SV0_QSO_tile68002_night20200315_1_merged.csv'
who = 'TMD'
VI_base    = 'desi-vi_SV0_QSO_tile68002_night20200315_1_'+who  # We add .csv to read in this file, and e.g. _z.png to print plots for this file
output_dir = VI_dir+'VI_check/'  # The directory in which we'll put the output. 
if not os.path.exists(output_dir):
  os.makedirs(output_dir)

m = pd.read_csv(VI_dir+VI_merged) 
v = pd.read_csv(VI_dir+VI_base+'.csv',delimiter = " , ", engine='python') #NOTE, prospect outputs data with spaces between commas, so you need the weird delimeter for prospect output, but not for normal csv output like we have done for the merged files
#print(m.keys())
#print(v.keys())

# If the VI agreed with Redrock, insert the Redrock results into the VI column.
v.loc[v['VI z']=='--', 'VI z'] = v.loc[v['VI z']=='--', 'Redrock z']
v['VI z'] = v['VI z'].astype(float)
v.loc[v['VI spectype']=='--', 'VI spectype'] = v.loc[v['VI spectype']=='--', 'Redrock spectype']

# Gather the results by TargetID, put them all in g, including only those items the VI'er looked at.
g = v.merge(m,how='left',on='TargetID')

# Find the disagreements
i_disagree_z = np.abs(g['VI z']-g['best z']) >=0.0033 
#print( g[['TargetID','best z','VI z']][i_disagree_z])
#print(sum(i_disagree_z))

i_disagree_class = abs(g['VI class']-g['best class']) >= 2  
#print( g[['TargetID','best class','VI class']][i_disagree_class])
#print(sum(i_disagree_class))

i_disagree_spectype = g['VI spectype'] != g['best spectype']
#print( g[['TargetID','best spectype','VI spectype']][i_disagree_spectype])
#print(sum(i_disagree_spectype))

i_good_bestz = g['best class'] >=2.5  # Choose only those redshifts that were given a quality >=2.5 in the merged results
i_good_VIz   = g['VI class']   >=3.0  # Choose only those redshifts that were given a quality >=3 by the VI
i_good_z = i_good_bestz | i_good_VIz

i_disagree_all = i_disagree_class | i_disagree_z | i_disagree_spectype
i_disagree_good = i_disagree_z & i_good_z
print('\nOut of %s objects, you disagreed with the merged results on %s redshifts, %s qualities, and %s spectypes.'%(len(g['TargetID']),sum(i_disagree_z),sum(i_disagree_class),sum(i_disagree_spectype)))
print('\nNumber of redshift disagreements on high-quality redshifts = %s.'%sum(i_disagree_good))
print('\nTotal number of disagreements = %s.'%sum(i_disagree_all))
print( g[['TargetID','best z','VI z','best class','VI class','best spectype','VI spectype']][i_disagree_all])

conflicts = ''
for target in g['TargetID'][i_disagree_all]:
  conflicts = conflicts+', '+str(target)
print('\nLook at your disagreements by opening the Prospect notebook viewer https://github.com/desihub/prospect/blob/master/doc/nb/Prospect_targetid.ipynb/ and cutting and pasting the following target list:')
print(conflicts[2:])

# Plot VI z vs Merged z and colour by Quality
plot_scatter_colors(g['best z'],g['VI z'],g['best class'],[0.0,4.0],'Best z','VI z','Quality (4 is good)',output_dir+VI_base+'_z.png',iextra=i_disagree_good,annotation='# of z disagreements on high-quality z: %s'%sum(i_disagree_good))
plot_histogram(g['best class'],np.arange(11)/2.-0.25,'Quality','Number of cases',output_dir+VI_base+'_quality.png',comp=g['VI class'],annotation='Blue=Merged VI; Orange=Single VI')

