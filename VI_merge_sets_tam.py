#!/usr/bin/env python
# coding: utf-8
# ## Identifying and resolving conflicts in a batch of VI files
# ## Adapted from Rita Tojero's jupyter notebook VI_merge_sets.ipynb
import os, sys, glob
import numpy as np
import inspect
from astropy.io import fits
from astropy.table import Table, join, vstack
from astropy.io import fits
import pandas as pd
import fnmatch
import matplotlib.pyplot as plt 
import re
from datetime import datetime
from pytz import timezone

# MAIN USER DEFINED INPUTS BELOW ARE "tiles", "nights", and "subset"; and importantly "output_name".
def display(x):
  print(x)

# Here you need to choose the tiles on which your objects were observed
tiledir   = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/'
tiles = ['68002']
nights = ['20200315']
petals = ['0','1', '2', '3', '4', '5', '6' ,'7', '8', '9']
subset = "_2_"  # YOU WANT TO CHANGE THIS EACH TIME, it defines "pattern" below
output_name = "desi-vi_SV0_QSO_tile"+tiles[0]+"_night"+nights[0]+subset+"merged"

on_nersc = True
if on_nersc:
  import desispec.io
  import desispec
  sys.path.append("/global/homes/t/tamarad/desihub/prospect/py")
  from prospect import utils_specviewer,plotframes

# Here you need to choose the tiles on which your objects were observed
  tiledir   = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/'
  tiles = ['68002']
  nights = ['20200315']
  petals = ['0','1', '2', '3', '4', '5', '6' ,'7', '8', '9']
  obs_db = utils_specviewer.make_targetdict(tiledir, petals=petals, tiles=tiles) # tiles = optional 

# Set to directory with all the VI files to merge
if on_nersc:
  VI_dir = os.environ['HOME']+'/SV/VI_files/SV0/QSO/'
else:
  #VI_dir = '/Users/uqtdavi1/Documents/Astro/DESI/SV_redshifts/SV0/QSO/'
  VI_dir = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/QSO/'


# We will read all the *.csv files in this directory. Change as needed.
all_files = os.listdir(VI_dir)
vi_files=[]

# Choose a subset
pattern = "desi*"+subset+"*.csv"
for entry in all_files:
  if fnmatch.fnmatch(entry, pattern):
    vi_files.append(entry)

# Prep the output files
output_file = VI_dir+'output/'+output_name+'.csv'
log_file = VI_dir+'output/'+output_name+'.log'
print(log_file)
# Make sure the path exists to write the output file
if not os.path.exists(VI_dir+'output'):
  os.makedirs(VI_dir+'output')
# Check whether a log file exists, if so copy it to a timestamped version in a subdirectory so we don't overwrite it
if os.path.isfile(log_file):
  print('Copying existing log file to subdirectory.')
  if not os.path.exists(VI_dir+'output/logs_old'):
    os.makedirs(VI_dir+'output/logs_old')
  os.system('cp '+log_file+' '+VI_dir+'output/logs_old/'+output_name+'.log_'+str(datetime.timestamp(datetime.now(timezone('UTC')))))

# Read the first file in to vi to set up vi
print(vi_files[0])
vi = pd.read_csv(VI_dir + vi_files[0], delimiter = " , ", engine='python')

# Read in the rest of the files and append them to vi
for i in range(1,len(vi_files)):
    print(vi_files[i])
    vi2 = pd.read_csv(VI_dir + vi_files[i], delimiter = " , ", engine='python')
    vi = vi.append(vi2, ignore_index=True)
    
# Change the column name to TARGETID to match standards elsewhere in DESI.
vi = vi.rename(columns={"TargetID": "TARGETID"})

#make groups of visual inspections, grouped by unique objects, and state number of single and multiple VIs
vi_gp = vi.groupby(['TARGETID'])
print('There are ' + str(len(vi)) + ' visual inspections of a total of ' + str(len(vi_gp)) + ' unique objects')

#vi is a dataframe
#pd.set_option('display.max_rows', 10) # This is for a notebook only
#display(vi.sort_values(by=['TARGETID']))
#vi.keys()

#----------------------------------------------------------------
# ### Merge with zbest files
# Add: fiberID, delta_chi2, flux information,.. anything else?
#read in fibermap info, loop over the files for all the petals
if on_nersc:
  tf = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[0])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='FIBERMAP')
  tspec = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[0])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='ZBEST')
  for i in range(1,len(petals)):
      tn = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[i])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='ZBEST')
      tnf = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[i])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='FIBERMAP')
      tspec = vstack([tspec,tn])
      tf = vstack([tf,tnf])

  tspec_df = tspec['TARGETID','DELTACHI2','ZWARN','ZERR'].to_pandas()
  tf_df = tf['TARGETID','FIBER','FLUX_G','FLUX_R','FLUX_Z','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z'].to_pandas()

  vi = vi.merge(tf_df, how='left', on='TARGETID')
  vi = vi.merge(tspec_df, how='left', on='TARGETID')


#----------------------------------------------------------------
# ### Adding a bunch of useful columns
#make new column with best redshift estimate for each VI - take VI redshift if available, else take Redrock redshift. 
#I am always assuming that the VI redshift, if provided, trumps over the Redrock redshift. 
vi['best z'] = vi['VI z']
vi.loc[vi['best z']=='--', 'best z'] = vi.loc[vi['best z']=='--', 'Redrock z']
vi['best z'] = vi['best z'].astype(float)

#add new column to find how much deviation there is in the redshift estimates 
vi['dz'] = vi.groupby('TARGETID')['best z'].transform(lambda x: ( (x.max() - x.min()) / (1+x.min()) ))

#if the deviation is small, fill best redshift with the mean redshift of the VI z (which may be a VI and a Redrock mean if only one VI changed the z)
vi.loc[vi['dz']< 0.0033,['best z']]=vi.loc[vi['dz']< 0.0033].groupby('TARGETID')['best z'].transform('mean')

#if the deviation is large, fill the best redshift with 999 so the merger has to enter a redshift
#vi.loc[vi['dz']>=0.0033,['best z']]=999.  # Turned off for now, so it can be used as a default.
#display(vi[['TARGETID','Redrock z','VI z','best z']].sort_values(by=['TARGETID']))

#make new column with best spectype estimate for each VI - take VI spectype if available, else take Redrock spectype 
#I am always assuming that the VI redshift, if provided, trumps over the Redrock redshift. 
vi['best spectype'] = vi['VI spectype']
vi.loc[vi['best spectype']=='--', 'best spectype'] = vi.loc[vi['best spectype']=='--', 'Redrock spectype']
#display(vi[['TARGETID','Redrock spectype','VI spectype','best spectype']].sort_values(by=['TARGETID']))

# Tam's edits
##make new column with issue flags - concatenate all issue flags from any VI.  We don't check these, we only check these if the VIs conflict for some other reason.
vi['all VI issues'] = vi.groupby('TARGETID')['VI issue'].transform(lambda x: ''.join(set(list(x))))
vi.loc[vi['all VI issues']!='--', 'all VI issues'] = vi.loc[vi['all VI issues']!='--', 'all VI issues'].transform(lambda x: ''.join(set(list(x))-set('-')))
#display(vi[['TARGETID','VI issue','all VI issues']].sort_values(by=['TARGETID']))

#add new columns, holding the mean of the flags and the maximum difference in flag classification
vi['best quality'] = vi.groupby('TARGETID')['VI class'].transform('mean')
# If the classification difference is too big, give it 999 so it is clear the merger has to assess it
vi['vi_class_diff'] = vi.groupby('TARGETID')['VI class'].transform(lambda x: (x.max()-x.min()))
vi.loc[vi['vi_class_diff']>1,['best quality']]=999
#display(vi[['TARGETID','VI class', 'vi_class_diff','best quality']].sort_values(by=['TARGETID']))

#add new column, with all comments concatenated
#vi['all VI comments'] = vi.groupby('TARGETID')['VI comment'].transform(lambda x: '|'.join(x))
vi['all VI comments'] = vi.groupby('TARGETID')['VI comment'].transform(lambda x: ' '.join(set(list(x))))
vi.loc[vi['all VI comments']!='--', 'all VI comments'] = vi.loc[vi['all VI comments']!='--', 'all VI comments'].transform(lambda x: x.replace("--",""))

#add new column, with the number of VI inspections for each object
vi['N_VI'] = vi.groupby('TARGETID')['TARGETID'].transform('count')

#add new column to hold comments from merger if needed
vi['merger comment'] = 'none'

#check all the new columns (keys) have been added correctly
#display(vi.sort_values(by=['TARGETID']))
#print(vi.keys())

#----------------------------------------------------------------
# Get a table that holds only the objects that have been inspected more than once, and for which the individual VI classifications differ by 2 or more, or delta z / (1 + z) > 0.0033, or there is disagreement in best spectype (these are the conflicts to resolve)
vi_gp = vi.groupby(['TARGETID'])
vi_conflict = vi_gp.filter(lambda x: ( ( (x['VI class'].max()-x['VI class'].min()) >= 2) 
                       | ( x['dz'].max()>=0.0033 ) 
                       | (not all(i == x['best spectype'].iloc[0] for i in x['best spectype'])) 
                       #| (not all(i == x['all VI issues'].iloc[0] for i in x['all VI issues']))  
                                     )
                       & (len(x) >= 2)) #x is a group by TARGETID

vi_badspectype = vi_gp.filter(lambda x: ( not (all(i == x['best spectype'].iloc[0] for i in x['best spectype'])) )
                       & (len(x) >= 2)) #x is a group by TARGETIDvi_badspectype


# Get the target IDs of the problematic objects and display in table form for a quick summary:
unique_targets = np.unique(vi_conflict['TARGETID'].tolist())
#print('Targets with problematic VI: ', unique_targets)
unique_target_csv = str(unique_targets[0])
for target in unique_targets[1:]:
  unique_target_csv = unique_target_csv+', '+str(target)
print('Copy the following list of problematic targets in to the "targets" list in Prospect_targetid.ipynb')
# On the wiki start from Computing/JupyterAtNERSC 
print('Targets with problematic VI: ', unique_target_csv)
print('Total number of conflicts to resolve: ', len(unique_targets))

#----------------------------------------------------------------
# ## This is where I resolve things manually - with care!!
# ### I think it's better to keep it in a notebook, as typos can be backtracked rather than a single manual edit of a text file
# 
# We edit either 'best quality', 'best redshift', 'best spectype, or 'all VI issues' to resolve conflict. At the end, we look for conflicts again and we should find none.
# 

#function to display the conflict in table format and open a prospect window
def display_conflict(conflict_id):
    #first, remind myself of the problem:
    display(vi[vi.TARGETID==unique_targets[conflict_id]][['TARGETID','Redrock spectype','VI spectype','best spectype','Redrock z','VI z','best z','VI class','best quality','VI issue','all VI issues','VI comment','merger comment','VI scanner']])
    #spectra, zcat= utils_specviewer.load_spectra_zcat_from_targets([unique_targets[conflict_id]], tiledir, obs_db)
    # VI interface in notebook
    #plotframes.plotspectra(spectra, zcatalog=zcat, title='Target_select', notebook=True, mask_type='CMX_TARGET',with_vi_widgets=False)

#first, keep a safe copy of the original dataframe
vi_safe = vi.copy()

# We will inspect each conflict on a prospect window, and resolve each conflict in turn
#copy this text to a new cell to display the conflict

def choose_spectype(argument):
  switcher = {
    's': 'STAR',
    'g': 'GALAXY',
    'q': 'QSO'
  }
  return switcher.get(argument,"Invalid_switch")


def issue_match(strg, search=re.compile(r'[^RCS-]').search):
  return not bool(search(strg))

print('-----------------------------------------------------------------')
print('Ready to do some merging?  Default values are in square brackets.')
print('-----------------------------------------------------------------')

# Start the counter
i=0
uselog='n'

# If a log file exists, read it in:
if os.path.isfile(log_file):
  log_old= pd.read_csv(log_file, delimiter = ', ', engine='python', comment='#')
  log_old = log_old.rename(columns={"TargetID": "TARGETID"}) # Temporary for reading in the old log files
  print('Log file:')
  display(log_old)
  nlog=len(log_old['index'])
  uselog=str(input("Log file exists with %s entries.  Press n to ignore, any other key to read in and continue:"%str(nlog)) or 'y')

# Prepare a new log file (this will overwrite the previous one, but a time-stamped copy of the previous one is saved in logs_old directory:
log=open(log_file,'w')
log.write('index, TARGETID, bestzmerge, bestclassmerge, bestspectypemerge, bestissuemerge, mergercomment\n')

if uselog != 'n':
  while i<nlog: 
    if unique_targets[i] != log_old.loc[i]['TARGETID']:
      print('WARNING: Log file and current input do not match.  Matching to TARGETID anyway.')
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'best z']         = log_old.loc[i]['bestzmerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'best quality']     = log_old.loc[i]['bestclassmerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'best spectype']  = log_old.loc[i]['bestspectypemerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'all VI issues']     = log_old.loc[i]['bestissuemerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'merger comment'] = log_old.loc[i]['mergercomment']
    log_text = str(i)+', '+str(log_old.loc[i]['TARGETID'])+', '+str(log_old.loc[i]['bestzmerge'])+', '+str(log_old.loc[i]['bestclassmerge'])+', '+log_old.loc[i]['bestspectypemerge']+', '+log_old.loc[i]['bestissuemerge']+', '+log_old.loc[i]['mergercomment']+'\n'
    log.write(log_text)
    i=i+1
  print('Read in log file, with %s entries, continuing from there.'%str(nlog))



while i<len(unique_targets): 
  print('test')
  print("%s/%s"%(i,len(unique_targets)-1))
  conflict = vi.loc[vi.TARGETID==unique_targets[i]]
  display_conflict(i)
  
  # Fix redshift
  tmp_bestz = conflict.loc[conflict.first_valid_index()]['best z']
  bestzmerge = float(input("Best z [%s]:"%str(tmp_bestz)) or tmp_bestz)
  
  # Fix class (quality)  #IMPROVEMENT: Make it not crash when someone enters a character
  tmp_bestclass = conflict.loc[conflict.first_valid_index()]['best quality']
  bestclassmerge = float(input("Best quality [%s]:"%str(tmp_bestclass)) or tmp_bestclass)
  while not(0<=bestclassmerge and 4>=bestclassmerge):
    print("Invalid choice. Quality must be between 0 and 4.")
    bestclassmerge = float(input("Best quality [%s]:"%str(tmp_bestclass)) or tmp_bestclass)

  # Fix spectype
  tmp_bestspectype = conflict.loc[conflict.first_valid_index()]['best spectype']
  sgq = str(input("Best spectype (s,g,q) [%s]:"%str(tmp_bestspectype)) or 'o')
  testchoice = (sgq=='s' or sgq=='g' or sgq=='q' or sgq=='o')
  while not(testchoice):
    print("Invalid choice, please choose s for STAR, g for GALAXY, q for QSO, or press enter to accept the default.")
    sgq = str(input("Best spectype (s,g,q) [%s]:"%str(tmp_bestspectype)) or 'o') 
    testchoice = (sgq=='s' or sgq=='g' or sgq=='q' or sgq=='o')
  if sgq == 'o':
    bestspectypemerge = tmp_bestspectype
  else:
    bestspectypemerge = choose_spectype(sgq)

  # Fix issue
  tmp_bestissue = conflict.loc[conflict.first_valid_index()]['all VI issues']
  bestissue = str(input("all VI issues (R=bad z, C=bad spectype, S=bad spectrum, --=no issue) [%s]:"%str(tmp_bestissue)) or 'o')
  testissue = (issue_match(bestissue) or bestissue=='o')
  while not(testissue):
    print("Invalid choice, please choose R for bad redshift, S for bad spectype, C for bad spectrum, or press enter to accept the default.")
    bestissue = str(input("all VI issues (R=bad z, S=bad spectype, C=bad spectrum) [%s]:"%str(tmp_bestissue)) or 'o')
    testissue = (issue_match(bestissue) or bestissue=='o')
  if bestissue == 'o':
    bestissuemerge = tmp_bestissue
  else:
    bestissuemerge = bestissue

  # Add comment
  mergercomment = str(input("Merger comment:") or 'None')

  print("Your results:  %s, %s, %s, %s, %s"%(str(bestzmerge), str(bestclassmerge), bestspectypemerge, bestissuemerge, mergercomment))
  happy = str(input("If you are NOT happy with these choices press n to enter them again, otherwise press any key to continue."))
  if happy != 'n':
    vi.loc[vi.TARGETID==unique_targets[i], 'best z'] = bestzmerge
    vi.loc[vi.TARGETID==unique_targets[i], 'best quality'] = bestclassmerge
    vi.loc[vi.TARGETID==unique_targets[i], 'best spectype'] = bestspectypemerge
    vi.loc[vi.TARGETID==unique_targets[i], 'all VI issues'] = bestissuemerge
    vi.loc[vi.TARGETID==unique_targets[i], 'merger comment'] = mergercomment

    log_text = str(i)+', '+str(unique_targets[i])+', '+str(bestzmerge)+', '+str(bestclassmerge)+', '+bestspectypemerge+', '+bestissuemerge+', '+mergercomment+'\n'
    print(log_text)
    log.write(log_text)
    print('----------------------------------------------------------------------')
    i=i+1

vi_unique = vi_gp['Redrock z', 'best z', 'best quality', 'Redrock spectype', 'best spectype', 'all VI issues', 'all VI comments', 'merger comment', 'N_VI'].first()
#print(vi_unique[vi_unique['best quality']<2.5]['best quality'].transform('count'))
print((vi_unique.loc[vi_unique['best quality']>=2.5,'best quality']).count())
print((vi_unique['best quality']).count())

number_of_total_objects = (vi_unique['best quality']).count()
number_of_good_redshifts = (vi_unique.loc[vi_unique['best quality']>=2.5,'best quality']).count()

print('Completed %s/%s checks out of %s objects (fraction needed checking=%s).'%(i,len(unique_targets),number_of_total_objects,len(unique_targets)/number_of_total_objects))
print('Resulting in %s/%s strong (>=2.5) redshifts for a success rate of %s.'%(number_of_good_redshifts,number_of_total_objects,number_of_good_redshifts/number_of_total_objects))
log.write('# Completed %s/%s checks out of %s objects (fraction needed checking=%s).\n'%(i,len(unique_targets),number_of_total_objects,len(unique_targets)/number_of_total_objects))
log.write('# Resulting in %s/%s strong (>=2.5) redshifts for a success rate of %s.\n'%(number_of_good_redshifts,number_of_total_objects,number_of_good_redshifts/number_of_total_objects))
log.close()

# ### and so on...
# 
# We should now recompute the conflicts, and not find any.

#vi_conflict = vi_gp.filter(lambda x: ( ( (x['VI class'].max()-x['VI class'].min()) >= 2) 
#                       | ( (x['best redshift'].max() - x['best redshift'].min()) / (1+x['best redshift'].min()) > 0.0033 ) 
#                       | (not all(i == x['best spectype'].iloc[0] for i in x['best spectype'])) )
#                       & (len(x) >= 2)) #x is a group by TARGETID
'''
vi_conflict_test = vi_gp.filter(lambda x: (all(i == 999 for i in x['best z'])) 
                               | (all(i==999 for i in x['best quality'])) 
                               & (len(x) >= 2) )
# This doesn't catch best spectype errors
pd.set_option('display.max_rows', 20)
display(vi_conflict_test[['TARGETID','Redrock z','VI z','best z','dz','best spectype','all VI issues','best quality','vi_class_diff','all VI comments']].sort_values(by=['TARGETID']))

unique_targets = np.unique(vi_conflict_test['TARGETID'].tolist())
print('Targets with problematic VI: ', unique_targets)
print('Total number of conflicts to resolve: ', len(unique_targets))


# Display anything that was missed (if "Total number of conflicts" isn't zero) and resolve!
for i in range(len(unique_targets)): 
    display(vi[vi.TARGETID==unique_targets[i]])


vi_conflict = vi_gp.filter(lambda x: ( ( (x['VI class'].max()-x['VI class'].min()) >= 2) 
                       | ( (x['best redshift'].max() - x['best redshift'].min()) / (1+x['best redshift'].min()) > 0.0033 ) 
                       | (not all(i == x['best spectype'].iloc[0] for i in x['best spectype'])) )
                       & (len(x) >= 2)) #x is a group by TARGETID


unique_targets = np.unique(vi_conflict['TARGETID'].tolist())
print('Targets with problematic VI: ', unique_targets)
print('Total number of conflicts to resolve: ', len(unique_targets))
'''

# ## Woohoo!

# ## Now we prepare to write to file. 
# 
# ### The important columns for the truth table construction are **best z**, **best quality**, **best spectype**, and **all VI issues**. 
print('Output to:',output_file)
print('Log to:',log_file)
if on_nersc:
#OLD  vi_gp['Redrock z', 'best z', 'best quality', 'Redrock spectype', 'best spectype', 'all VI issues', 'all VI comments', 'merger comment', 'N_VI', 'DELTACHI2','FIBER','FLUX_G','FLUX_R','FLUX_Z','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z'].first().to_csv(output_file)
  vi_gp['Redrock z', 'best z', 'best quality', 'Redrock spectype', 'best spectype', 'all VI issues', 'all VI comments', 'merger comment','N_VI','DELTACHI2', 'ZWARN', 'ZERR','FIBER','FLUX_G', 'FLUX_R', 'FLUX_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MW_TRANSMISSION_G','MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z'].first().to_csv(output_file)
else:
  vi_gp['Redrock z', 'best z', 'best quality', 'Redrock spectype', 'best spectype', 'all VI issues', 'all VI comments', 'merger comment', 'N_VI'].first().to_csv(output_file)


# Check that merged file reads in OK - check comments
#merged_file = pd.read_csv(output_file)
#merged_file.keys()
#merged_file
#merged_file.loc[9]['merger comment']
#merged_file.loc[9]['all VI comments']

