#!/usr/bin/env python
# coding: utf-8
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

VI_dir = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO/merged/'
tiles = ['68002']
nights = ['20200315']
#subset = "_3_"

input_file       = VI_dir+"desi-vi_SV0_QSO_tile"+tiles[0]+"_night"+nights[0]+"_merged_all.csv"
output_file_base = VI_dir+"desi-vi_SV0_QSO_"+tiles[0]+"_"+nights[0]
m = pd.read_csv(input_file)
print(m.keys())

def plot_scatter_colors(x,y,col,clim,xlabel,ylabel,collabel,plot_name,xlim=0,ylim=0,psize=8.0,annotation=''):
  print('Plotting scatter with colours')
  plt.figure(dpi=200)
  plt.scatter(x,y, s=psize,
              c=col, cmap=plt.cm.get_cmap('jet_r'))
  plt.colorbar(ticks=range(5), label=collabel)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if xlim !=0: plt.xlim(xlim)
  if ylim !=0: plt.ylim(ylim)
  plt.annotate(annotation,xy=(0.1,0.1),xycoords='axes fraction')
  plt.clim(clim[0],clim[1])
  plt.savefig(plot_name,bbox_inches='tight')
  plt.close()

def plot_histogram(x,bin_edges,xlabel,ylabel,annotation,title,plot_name):
  plt.figure(dpi=200)
  plt.hist(x, bins = bin_edges,rwidth=0.85)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  plt.annotate(annotation,xy=(0.1,0.9),xycoords='axes fraction')
  plt.savefig(plot_name,bbox_inches='tight')
  plt.close()

# Enumerate the spectypes so they can be easily plotted
ispectype = np.zeros(len(m['best spectype']))
ispectype[m['best spectype']=='STAR'] = int(0)
ispectype[m['best spectype']=='GALAXY'] = int(1)
ispectype[m['best spectype']=='QSO'] = int(2)
print('len(ispectype)=',len(ispectype))
print(len(ispectype[m['best spectype']=='STAR']))
print(len(ispectype[m['best spectype']=='GALAXY']))
print(len(ispectype[m['best spectype']=='QSO']))

print('S' in m['best issue'])
exit()

# Plot Best z vs Redrock z and colour by Quality
plot_scatter_colors(m['Redrock z'],m['best z'],m['best class'],[0.0,4.0],'Redrock z','Best z','Quality (4 is good)',output_file_base+'_bestz-vs-Redrock.png')

# Plot Best Class (quality) vs Flux and colour by Spectype
plot_scatter_colors(m['FLUX_R'],m['best class']+(np.random.randn(len(m['best class']))-1.0)*0.005,ispectype,[0.0,2.0],
  'FLUX_R','Quality (4 is good)','Spectype (0=STAR,1=GALAXY,2=QSO)',output_file_base+'_quality-vs-flux.png',xlim=[0,20],psize=4.0,annotation='Slight random added to Quality for visibility')

# Calculate and print number of good redshifts
ngood = np.sum((m['best class']>=2.5))
nbad  = len(m['best class'])-np.sum((m['best class']>=2.5))
percentfail = float(nbad)/(float(ngood)+float(nbad))*100.
print('Good redshifts=',ngood)
print('Bad redshifts =',nbad)
print('Fail % = ',percentfail)
#25 have q=2.5
print(25./(float(ngood)+float(nbad)))

# Plot histogram of classifications (qualities)
plot_histogram(m['best class'],[-0.5,0.5,1.5,2.5,3.5,4.5],'Quality (4 is good)','Number of VI results',
  'Percentage failure (<2.5) = {:.1f}%'.format(percentfail),
  'All VI results QSO SV0 tile='+tiles[0]+' night='+nights[0],output_file_base+'_histogram.png')

# Select only those classified as good.
good = m[m['best class']>=2.5]
ispectypegood = ispectype[m['best class']>=2.5]
#print(good['best class'])
bad = m[m['best class']<2.5]

print('NSTAR=',len(ispectypegood[ispectypegood==0]))
print('NGALAXY=',len(ispectypegood[ispectypegood==1]))
print('NQSO=',len(ispectypegood[ispectypegood==2]))


# Plot Best z vs Redrock z and colour by Quality for good redshifts only
plot_scatter_colors(good['Redrock z'],good['best z'],good['best class'],[0.0,4.0],'Redrock z','Best z','Quality (4 is good)',output_file_base+'_bestz-vs-Redrock_good.png')
# Plot Best z vs Redrock z and colour by Quality for bad redshifts only
plot_scatter_colors(bad['Redrock z'],bad['best z'],bad['best class'],[0.0,4.0],'Redrock z','Best z','Quality (4 is good)',output_file_base+'_bestz-vs-Redrock_bad.png')

# Find good redshifts where Redrock disagreed with best VI
good_failures = good[ np.abs(good['Redrock z']-good['best z'])/good['best z']>0.0033 ] 
#print(good_failures[['best z','Redrock z']])

# Plot histogram of good failures
plot_histogram(good_failures['best z'],[-0.05,0.05,0.5,1.0,1.5, 2.0,2.5, 3.0,3.5,4.0, 4.5],'Redshift','Number of VI results',
  'Secure Redshifts Redrock misidentified',
  'All VI results QSO SV0 tile='+tiles[0]+' night='+nights[0],output_file_base+'_histogram_redshifts_good_failures.png')
print(len(good_failures),len(good),len(good_failures)/len(good))

# Plot histogram of fluxes of successful redshifts
flux_bins = np.arange(31)
plot_histogram(good['FLUX_R'],flux_bins,'FLUX_R','Number successful redshifts','Total number of successful redshifts = %s'%ngood,
  'All VI results QSO SV0 tile='+tiles[0]+' night='+nights[0],output_file_base+'_histogram_FLUXR_success.png')
plot_histogram(bad['FLUX_R'],flux_bins,'FLUX_R','Number UNsuccessful redshifts','Total number of unsuccessful redshifts = %s'%nbad,
  'All VI results QSO SV0 tile='+tiles[0]+' night='+nights[0],output_file_base+'_histogram_FLUXR_unsuccess.png')

# Plot fraction of successful redshifts vs flux
flux_bins = np.arange(21)
goodhist,bin_edges = np.histogram(good['FLUX_R'],bins=flux_bins)
badhist,bin_edges = np.histogram(bad['FLUX_R'],bins=flux_bins)
totalhist,bin_edges = np.histogram(m['FLUX_R'],bins=flux_bins)
totalhist[totalhist==0] = 1.e5

goodfrac = goodhist/totalhist
badfrac = badhist/totalhist
bin_centres=(bin_edges[0:-1]+bin_edges[1:])/2
plt.figure(dpi=200)
plt.scatter(bin_centres,badfrac,marker='x',c='red')#,color='red',label='failure (quality<2.5)')
plt.scatter(bin_centres,goodfrac,marker='o',c='green')#,'o',color='green',label='success (quality>=2.5)')
plt.annotate('Successful   (quality>=2.5)',xy=(0.5,0.9),xycoords='axes fraction',color='green')
plt.annotate('Unsuccessful (quality< 2.5)',xy=(0.5,0.1),xycoords='axes fraction',color='red')
plt.xlabel('FLUX_R')
plt.ylabel('Fraction success (o) or failure (x)')
plt.savefig(output_file_base+'_FracSuccess_vs_FLUXR.png',bbox_inches='tight')
plt.close()



# Plot Best z vs Redrock z and colour by Quality for good redshifts only
plot_scatter_colors(m['FLUX_R'],m['FLUX_G'],m['best class'],[0.0,4.0],'FLUX_R','FLUX_G','Quality (4 is good)',output_file_base+'_FLUXR_vs_FLUXG.png',xlim=[0,20],ylim=[0,20])
plot_scatter_colors(m['FLUX_R'],m['FLUX_G'],m['best class'],[0.0,4.0],'FLUX_R','FLUX_G','Quality (4 is good)',output_file_base+'_FLUXR_vs_FLUXG_zoom.png',xlim=[0,4],ylim=[0,4])
plot_scatter_colors(m['FLUX_R'],m['FLUX_Z'],m['best class'],[0.0,4.0],'FLUX_R','FLUX_Z','Quality (4 is good)',output_file_base+'_FLUXR_vs_FLUXZ.png',xlim=[0,20],ylim=[0,20])
plot_scatter_colors(m['FLUX_R'],m['FLUX_Z'],m['best class'],[0.0,4.0],'FLUX_R','FLUX_Z','Quality (4 is good)',output_file_base+'_FLUXR_vs_FLUXZ_zoom.png',xlim=[0,4],ylim=[0,6])
plot_scatter_colors(m['FLUX_R'],m['FIBERFLUX_R'],m['best class'],[0.0,4.0],'FLUX_R','FIBERFLUX_R','Quality (4 is good)',output_file_base+'_FLUXR_vs_FIBERFLUXR.png',xlim=[0,20],ylim=[0,15])
plot_scatter_colors(m['FLUX_Z'],m['FIBERFLUX_Z'],m['best class'],[0.0,4.0],'FLUX_Z','FIBERFLUX_Z','Quality (4 is good)',output_file_base+'_FLUXZ_vs_FIBERFLUXZ.png',xlim=[0,20],ylim=[0,15])
plot_scatter_colors(m['FLUX_G'],m['FIBERFLUX_G'],m['best class'],[0.0,4.0],'FLUX_G','FIBERFLUX_G','Quality (4 is good)',output_file_base+'_FLUXG_vs_FIBERFLUXG.png',xlim=[0,20],ylim=[0,15])

plot_scatter_colors(good['FLUX_R'],good['best z'],good['best class'],[0.0,4.0],'FLUX_R','Best z','Quality (4 is good)',output_file_base+'_FLUXR_vs_Bestz.png',xlim=[0,20])
plot_scatter_colors(good['FLUX_Z'],good['best z'],good['best class'],[0.0,4.0],'FLUX_Z','Best z','Quality (4 is good)',output_file_base+'_FLUXZ_vs_Bestz.png',xlim=[0,20])
plot_scatter_colors(good['FLUX_G'],good['best z'],good['best class'],[0.0,4.0],'FLUX_G','Best z','Quality (4 is good)',output_file_base+'_FLUXG_vs_Bestz.png',xlim=[0,20])

plot_scatter_colors(good['FLUX_R'],good['best z'],ispectypegood,[0.0,2.0],'FLUX_R','Best z','Spectype (0=STAR,1=GALAXY,2=QSO)',output_file_base+'_FLUXR_vs_Bestz_types.png',xlim=[0,30])
plot_scatter_colors(good['FLUX_Z'],good['best z'],ispectypegood,[0.0,2.0],'FLUX_Z','Best z','Spectype (0=STAR,1=GALAXY,2=QSO)',output_file_base+'_FLUXZ_vs_Bestz_types.png',xlim=[0,30])
plot_scatter_colors(good['FLUX_G'],good['best z'],ispectypegood,[0.0,2.0],'FLUX_G','Best z','Spectype (0=STAR,1=GALAXY,2=QSO)',output_file_base+'_FLUXG_vs_Bestz_types.png',xlim=[0,30])


exit()

# PLOT ALL DATA
if False:
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






