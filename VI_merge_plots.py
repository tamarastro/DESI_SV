#!/usr/bin/env python
# coding: utf-8
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

VI_dir = '/Users/uqtdavi1/Documents/Astro/DESI/SV_redshifts/SV0/merged/'
tiles = ['68002']
nights = ['20200315']
subset = "_3_"

input_file = VI_dir+"desi-vi_SV0_QSO_tile"+tiles[0]+"_night"+nights[0]+subset+"merged.csv"
print(input_file)

if False:
  m = pd.read_csv(input_file)
  print(m.keys())
  print(len(m['TargetID']))
  ax = plt.subplot(111)
    
  for i in np.arange(len(m['TargetID'])):
    print(m[['Redrock z','best z','best class']])

    #ax.plot(m['Redrock z'][i],m['best z'][i],'o',color=cmap)
    #fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    #view_colormap('viridis')

  plt.scatter(m['Redrock z'],m['best z'], s=16.0,
              c=m['best class'], cmap=plt.cm.get_cmap('jet_r'))
  plt.colorbar(ticks=range(5), label='Quality (4 is good)')
  plt.xlabel('Redrock z')
  plt.ylabel('Best z')
  plt.clim(0.0, 4.0)

  output_file =input_file[0:-4]+'.png'
  plt.savefig(output_file,bbox_inches='tight')
  plt.close()

# PLOT ALL DATA
if True:
  input_file = VI_dir+"desi-vi_SV0_QSO_tile"+tiles[0]+"_night"+nights[0]+"all.csv"
  m = pd.read_csv(input_file)
  print(m.keys())
  print(len(m['TargetID']))
  #for i in np.arange(len(m['TargetID'])):
  #print(m[['Redrock z','VI z','best z','VI class']])
  m.loc[m['best z']==999,'best z']=m.loc[m['best z']==999,'Redrock z']
  #print(m[['Redrock z','VI z','best z','VI class']])

  plt.scatter(m['Redrock z'],m['best z'], s=16.0,
    c=m['VI class'], cmap=plt.cm.get_cmap('jet_r'))
  plt.colorbar(ticks=range(5), label='Quality (4 is good)')
  plt.xlabel('Redrock z')
  plt.ylabel('VI z')
  plt.clim(0.0, 4.0)
  output_file =input_file[0:-4]+'.png'
  plt.savefig(output_file,bbox_inches='tight')
  plt.close()

  plt.hist(m['VI class'], bins = [-0.5,0.5,1.5,2.5,3.5,4.5],rwidth=0.85)
  plt.xlabel('Quality (4 is good)')
  plt.ylabel('Number of VI results')
  plt.title('All VI results QSO SV0 tile='+tiles[0]+' night='+nights[0])
  ngood = np.sum((m['VI class']>2.5))
  nbad  = len(m['VI class'])-np.sum((m['VI class']>2.5))
  print('Good redshifts=',ngood)
  print('Bad redshifts =',nbad)
  print('Fail % = ',float(nbad)/float(ngood))
  output_file =input_file[0:-4]+'hist.png'
  plt.savefig(output_file,bbox_inches='tight')
  plt.close()


  plt.scatter(m['Redrock z'],m['VI class'], s=16.0,
    c=m['VI class'], cmap=plt.cm.get_cmap('jet_r'))
  plt.colorbar(ticks=range(5), label='Quality (4 is good)')
  plt.xlabel('Redrock z')
  plt.ylabel('VI class')
  plt.clim(0.0, 4.0)
  output_file =input_file[0:-4]+'qualityvsz.png'
  plt.savefig(output_file,bbox_inches='tight')
  plt.close()






