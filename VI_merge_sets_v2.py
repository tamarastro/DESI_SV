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
from VI_merge_functions_v2 import *
pd.options.display.max_colwidth = 150

on_nersc = False 
if on_nersc:
  import desispec.io
  import desispec

#--------------------------------------------------------------------------------------------------
# USER DEFINED INPUTS BELOW ARE "tiles", "nights", and "subset"; and importantly "output_name".
#--------------------------------------------------------------------------------------------------
if on_nersc:  
  # Set to directory with all the VI files to merge
  VI_dir = os.environ['HOME']+'/SV/VI_files/SV0/Blanc/QSO/'
else:
  VI_dir = '/Users/uqtdavi1/Documents/programs/DESI/SV/VI_files/SV0/Blanc/BGS/'

# Here you need to choose the tiles on which your objects were observed
tiledir   = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/'
tiles = ['80613']
nights = ['20201215']
petals = ['0','1', '2', '3', '4', '5', '6' ,'7', '8', '9']
subset = "_2_"  # YOU WANT TO CHANGE THIS EACH TIME, it defines "pattern" below.  Set to "" to use all.
output_name = "desi-vi_BGS_tile"+tiles[0]+"_nightdeep"+subset+"merged"

# Prep the output files
output_file = VI_dir+'output/'+output_name+'.csv'
log_file = VI_dir+'output/'+output_name+'.log'
#--------------------------------------------------------------------------------------------------
# END USER DEFINED INPUTS
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Read in all the data and combine the files ready for merging
#--------------------------------------------------------------------------------------------------
# Read in the data
vi = read_in_data(VI_dir,tiles[0],subset)

#make groups of visual inspections, grouped by unique objects, and state number of single and multiple VIs
vi_gp = vi.groupby(['TARGETID'])
print('There are ' + str(len(vi)) + ' visual inspections of a total of ' + str(len(vi_gp)) + ' unique objects')

#----------------------------------------------------------------
# Add extra data from zbest and fibermap files
# Add: fiberID, delta_chi2, flux information, MW transmission
# read in fibermap info, loop over the files for all the petals
if False: #on_nersc:  #Change back when we have the auxiliary data
  vi = add_auxiliary_data(vi,tiledir,tiles,nights,petals)

#----------------------------------------------------------------
# ### Adding useful columns.  Automatically decide on best_z, spectype, and quality where possible.  Concatenate issues and comments.
choose_best_z(vi)
choose_best_spectype(vi)
choose_best_quality(vi)
concatenate_all_issues(vi)
concatenate_all_comments(vi)
add_extra_details(vi)

#check all the new columns (keys) have been added correctly
print(vi.keys())

#----------------------------------------------------------------
# Get a table that holds only the objects that have been inspected more than once, and for which the individual VI classifications differ by 2 or more, or delta z / (1 + z) > 0.0033, or there is disagreement in best_spectype (these are the conflicts to resolve)
vi_gp = vi.groupby(['TARGETID'])
vi_conflict = find_conflicts(vi_gp)

# Get the target IDs of the problematic objects and display list that can be read in to Prospect_targetid.ipynb:
unique_targets = np.unique(vi_conflict['TARGETID'].tolist())
print_conflicts_for_prospect(unique_targets)


#----------------------------------------------------------------
# ## Resolve conflicts manually
#----------------------------------------------------------------
# We edit either 'best_quality', 'best redshift', 'best_spectype, or 'all_VI_issues' to resolve conflict. At the end, we look for conflicts again and we should find none.

print('-----------------------------------------------------------------')
print('Ready to do some merging?  Default values are in square brackets.')
print('-----------------------------------------------------------------')

# Start the counter
i=0
uselog='n'

#----------------------------------
# Check whether a log file exists, if so copy it to a timestamped version in a subdirectory so we don't overwrite it
if os.path.isfile(log_file):
  print('Copying existing log file to subdirectory.')
  if not os.path.exists(VI_dir+'output/logs_old'):
    os.makedirs(VI_dir+'output/logs_old')
  os.system('cp '+log_file+' '+VI_dir+'output/logs_old/'+output_name+'.log_'+str(datetime.timestamp(datetime.now(timezone('UTC')))))

# If a log file exists, read it in:
if os.path.isfile(log_file):
  log_old= pd.read_csv(log_file, delimiter = ', ', engine='python', comment='#', keep_default_na=False)
  log_old = log_old.rename(columns={"TargetID": "TARGETID"}) # Temporary for reading in the old log files
  print('Log file:')
  print(log_old)
  nlog=len(log_old['index'])
  uselog=str(input("Log file exists with %s entries.  Press n to ignore, any other key to read in and continue:"%str(nlog)) or 'y')


# Prepare a new log file (this will overwrite the previous one, but a time-stamped copy of the previous one is saved in logs_old directory:
log=open(log_file,'w')
log.write('index, TARGETID, bestzmerge, bestclassmerge, bestspectypemerge, bestissuemerge, mergercomment\n')

if uselog != 'n':
  while i<nlog: 
    if unique_targets[i] != log_old.loc[i]['TARGETID']:
      print('WARNING: Log file and current input do not match.  Matching to TARGETID anyway.')
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'best_z']         = log_old.loc[i]['bestzmerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'best_quality']   = log_old.loc[i]['bestclassmerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'best_spectype']  = log_old.loc[i]['bestspectypemerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'all_VI_issues']  = log_old.loc[i]['bestissuemerge']
    vi.loc[vi.TARGETID==log_old.loc[i]['TARGETID'], 'merger_comment'] = log_old.loc[i]['mergercomment']
    log_text = str(i)+', '+str(log_old.loc[i]['TARGETID'])+', '+str(log_old.loc[i]['bestzmerge'])+', '+str(log_old.loc[i]['bestclassmerge'])+', '+log_old.loc[i]['bestspectypemerge']+', '+log_old.loc[i]['bestissuemerge']+', '+log_old.loc[i]['mergercomment']+'\n'
    log.write(log_text)
    i=i+1
  print('Read in log file, with %s entries, continuing from there.'%str(nlog))
#----------------------------------

#----------------------------------
# Perform the merging if no log file was used, or continue from where we left off if a log file partially filled in the answers
while i<len(unique_targets): 
  print('test')
  print("%s/%s"%(i,len(unique_targets)-1))
  conflict = vi.loc[vi.TARGETID==unique_targets[i]]
  print_conflict(vi,unique_targets,i)
  
  # Fix redshift
  tmp_bestz = conflict.loc[conflict.first_valid_index()]['best_z']
  bestzmerge = float(input("best_z [%s]:"%str(tmp_bestz)) or tmp_bestz)
  
  # Fix class (quality)  #IMPROVEMENT NEEDED: Make it not crash when someone enters a character
  tmp_bestclass = conflict.loc[conflict.first_valid_index()]['best_quality']
  bestclassmerge = float(input("best_quality [%s]:"%str(tmp_bestclass)) or tmp_bestclass)
  while not(0<=bestclassmerge and 4>=bestclassmerge):
    print("Invalid choice. Quality must be between 0 and 4.")
    bestclassmerge = float(input("best_quality [%s]:"%str(tmp_bestclass)) or tmp_bestclass)

  # Fix spectype
  tmp_bestspectype = conflict.loc[conflict.first_valid_index()]['best_spectype']
  sgq = str(input("best_spectype (s,g,q) [%s]:"%str(tmp_bestspectype)) or 'o')
  testchoice = (sgq=='s' or sgq=='g' or sgq=='q' or sgq=='o')
  while not(testchoice):
    print("Invalid choice, please choose s for STAR, g for GALAXY, q for QSO, or press enter to accept the default.")
    sgq = str(input("best_spectype (s,g,q) [%s]:"%str(tmp_bestspectype)) or 'o') 
    testchoice = (sgq=='s' or sgq=='g' or sgq=='q' or sgq=='o')
  if sgq == 'o':
    bestspectypemerge = tmp_bestspectype
  else:
    bestspectypemerge = choose_spectype(sgq)

  # Fix issue
  tmp_bestissue = conflict.loc[conflict.first_valid_index()]['all_VI_issues']
  bestissue = str(input("all_VI_issues (R=bad z, C=bad spectype, S=bad spectrum, N=no issue) [%s]:"%str(tmp_bestissue)) or 'o')
  testissue = (issue_match(bestissue) or bestissue=='o')
  while not(testissue):
    print("Invalid choice, please choose R for bad redshift, S for bad spectype, C for bad spectrum, N to report no issue, or press enter to accept the default.")
    bestissue = str(input("all_VI_issues (R=bad z, S=bad spectype, C=bad spectrum) [%s]:"%str(tmp_bestissue)) or 'o')
    testissue = (issue_match(bestissue) or bestissue=='o')
  if bestissue == 'N':
    bestissue = ''
  if bestissue == 'o':
    bestissuemerge = tmp_bestissue
  else:
    bestissuemerge = bestissue

  # Add comment
  mergercomment = str(input("merger_comment:") or 'None')

  print("Your results:  %s, %s, %s, %s, %s"%(str(bestzmerge), str(bestclassmerge), bestspectypemerge, bestissuemerge, mergercomment))
  happy = str(input("If you are NOT happy with these choices press n to enter them again, otherwise press any key to continue."))
  if happy != 'n':
    vi.loc[vi.TARGETID==unique_targets[i], 'best_z'] = bestzmerge
    vi.loc[vi.TARGETID==unique_targets[i], 'best_quality'] = bestclassmerge
    vi.loc[vi.TARGETID==unique_targets[i], 'best_spectype'] = bestspectypemerge
    vi.loc[vi.TARGETID==unique_targets[i], 'all_VI_issues'] = bestissuemerge
    vi.loc[vi.TARGETID==unique_targets[i], 'merger_comment'] = mergercomment

    log_text = str(i)+', '+str(unique_targets[i])+', '+str(bestzmerge)+', '+str(bestclassmerge)+', '+bestspectypemerge+', '+bestissuemerge+', '+mergercomment+'\n'
    print(log_text)
    log.write(log_text)
    print('----------------------------------------------------------------------')
    i=i+1

#----------------------------------
# Now they've been given best values, make a new vi with only one instance of each TARGETID
vi_unique = vi_gp['Redrock_z', 'best_z', 'best_quality', 'Redrock_spectype', 'best_spectype', 'all_VI_issues', 'all_VI_comments', 'merger_comment', 'N_VI'].first()

number_of_total_objects = (vi_unique['best_quality']).count()
number_of_good_redshifts = (vi_unique.loc[vi_unique['best_quality']>=2.5,'best_quality']).count()

print('Completed %s/%s checks out of %s objects (fraction needed checking=%s).'%(i,len(unique_targets),number_of_total_objects,len(unique_targets)/number_of_total_objects))
print('Resulting in %s/%s strong (>=2.5) redshifts for a success rate of %s.'%(number_of_good_redshifts,number_of_total_objects,number_of_good_redshifts/number_of_total_objects))
log.write('# Completed %s/%s checks out of %s objects (fraction needed checking=%s).\n'%(i,len(unique_targets),number_of_total_objects,len(unique_targets)/number_of_total_objects))
log.write('# Resulting in %s/%s strong (>=2.5) redshifts for a success rate of %s.\n'%(number_of_good_redshifts,number_of_total_objects,number_of_good_redshifts/number_of_total_objects))
log.close()

#----------------------------------
# ## Write to file. 
# 
# ### The important columns for the truth table construction are **best_z**, **best_quality**, **best_spectype**, and **all_VI_issues**. 
print('\nOutput to:',output_file)
print('Log to:',log_file)
if False: #on_nersc:
  print_merged_file(vi_gp,output_file)
else:
  vi_gp['Redrock_z', 'best_z', 'best_quality', 'Redrock_spectype', 'best_spectype', 'all_VI_issues', 'all_VI_comments', 'merger_comment', 'N_VI'].first().to_csv(output_file)


