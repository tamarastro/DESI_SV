# VI_merge_functions_v1
# File for VI merging functions for both notebook and command-line use.
#
import os, sys, glob
import fnmatch
import re
from astropy.table import Table, join, vstack
import pandas as pd

def choose_best_z(vi):
  #make new column with best redshift estimate for each VI - take VI redshift if available, else take Redrock redshift. 
  vi['best z'] = vi['VI z']
  vi.loc[vi['best z']=='--', 'best z'] = vi.loc[vi['best z']=='--', 'Redrock z']
  vi['best z'] = vi['best z'].astype(float)

  #add new column to find how much deviation there is in the redshift estimates 
  vi['dz'] = vi.groupby('TARGETID')['best z'].transform(lambda x: ( (x.max() - x.min()) / (1+x.min()) ))

  #if the deviation is small, fill best redshift with the mean redshift of the VI z (which may be a VI and a Redrock mean if only one VI changed the z)
  vi.loc[vi['dz']< 0.0033,['best z']]=vi.loc[vi['dz']< 0.0033].groupby('TARGETID')['best z'].transform('mean')

def choose_best_spectype(vi):
  #make new column with best spectype estimate for each VI - take VI spectype if available, else take Redrock spectype 
  vi['best spectype'] = vi['VI spectype']
  vi.loc[vi['best spectype']=='--', 'best spectype'] = vi.loc[vi['best spectype']=='--', 'Redrock spectype']

def choose_best_quality(vi):
  #add new columns, holding the mean of the flags and the maximum difference in flag classification
  vi['best quality'] = vi.groupby('TARGETID')['VI class'].transform('mean')
  # If the classification difference is too big, give it 999 so it is clear the merger has to assess it
  vi['vi_class_diff'] = vi.groupby('TARGETID')['VI class'].transform(lambda x: (x.max()-x.min()))
  vi.loc[vi['vi_class_diff']>1,['best quality']]=999

def concatenate_all_issues(vi):
  ##make new column with issue flags - concatenate all issue flags from any VI.  We don't check these, we only check these if the VIs conflict for some other reason.
  vi['all VI issues'] = vi.groupby('TARGETID')['VI issue'].transform(lambda x: ''.join(set(list(x))))
  vi.loc[vi['all VI issues']!='--', 'all VI issues'] = vi.loc[vi['all VI issues']!='--', 'all VI issues'].transform(lambda x: ''.join(set(list(x))-set('-')))

def concatenate_all_comments(vi):
  #add new column, with all comments concatenated
  vi['all VI comments'] = vi.groupby('TARGETID')['VI comment'].transform(lambda x: ' '.join(set(list(x))))
  vi.loc[vi['all VI comments']!='--', 'all VI comments'] = vi.loc[vi['all VI comments']!='--', 'all VI comments'].transform(lambda x: x.replace("--",""))
  vi['all VI comments'] = vi['all VI comments'].transform(lambda x: x.strip())

def add_extra_details(vi):
  #add new column, with the number of VI inspections for each object
  vi['N_VI'] = vi.groupby('TARGETID')['TARGETID'].transform('count')

  #add new column to hold comments from merger if needed
  vi['merger comment'] = 'none'

def read_in_data(VI_dir,tile,subset):
  # We will read all the *.csv files in this directory. Change as needed.
  all_files = os.listdir(VI_dir)
  vi_files=[]

  # Choose a subset
  pattern = "desi*"+tile+"*"+subset+"*.csv"
  for entry in all_files:
    if fnmatch.fnmatch(entry, pattern):
      vi_files.append(entry)

  # Make sure the path exists to write the output file
  if not os.path.exists(VI_dir+'output'):
    os.makedirs(VI_dir+'output')

  # Read the first file in to vi to set up vi
  print('VI Files:')
  print(vi_files[0])
  vi = pd.read_csv(VI_dir + vi_files[0], delimiter = " , ", engine='python')

  # Read in the rest of the files and append them to vi
  for i in range(1,len(vi_files)):
      print(vi_files[i])
      vi2 = pd.read_csv(VI_dir + vi_files[i], delimiter = " , ", engine='python')
      vi = vi.append(vi2, ignore_index=True)
      
  # Change the column name to TARGETID to match standards elsewhere in DESI.
  #vi = vi.rename(columns={"TargetID": "TARGETID"})
  return vi

def add_auxiliary_data(vi,tiledir,tiles,nights,petals):
  tf = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[0])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='FIBERMAP')
  tspec = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[0])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='ZBEST')
  for i in range(1,len(petals)):
      tfn = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[i])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='FIBERMAP')
      tf = vstack([tf,tfn])
      tspecn = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[i])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='ZBEST')
      tspec = vstack([tspec,tspecn])
      
  tf_df = tf['TARGETID','FIBER','FLUX_G','FLUX_R','FLUX_Z','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z'].to_pandas()
  tspec_df = tspec['TARGETID','DELTACHI2','ZWARN','ZERR'].to_pandas()
  
  vi = vi.merge(tf_df, how='left', on='TARGETID')
  vi = vi.merge(tspec_df, how='left', on='TARGETID')
  return vi

def find_conflicts(vi_gp):
  # Choose spectra where VI has disagreed
  vi_conflict = vi_gp.filter(lambda x: (                      #x is a group by TARGETID
    ( ( x['VI class'].max()-x['VI class'].min()) >= 2     )   # Quality differs by 2 or more.
    | ( x['dz'].max()                            >=0.0033 )   # Redshift difference is >=0.0033 (approx 1000km/s at low-z).
    | ( not all(i == x['best spectype'].iloc[0] for i in x['best spectype']) )  # Best spectype differs.
    )
    & (len(x) >= 2)) # And there is more than one VI
  return(vi_conflict)

def print_conflict(conflict_id):
  #function to display the conflict in table format and open a prospect window
  print(vi[vi.TARGETID==unique_targets[conflict_id]][['TARGETID','Redrock spectype','VI spectype','best spectype','Redrock z','VI z','best z','VI class','best quality','VI issue','all VI issues','VI comment','merger comment','VI scanner']])

def choose_spectype(argument):
  switcher = {
    's': 'STAR',
    'g': 'GALAXY',
    'q': 'QSO'
  }
  return switcher.get(argument,"Invalid_switch")

def issue_match(strg, search=re.compile(r'[^RCS-]').search):
  # Clean up issues
  return not bool(search(strg))

def print_conflicts_for_prospect(unique_targets):
  unique_target_csv = str(unique_targets[0])
  for target in unique_targets[1:]:
    unique_target_csv = unique_target_csv+', '+str(target)
  print('Copy the following list of problematic targets in to the "targets" list in Prospect_targetid.ipynb')
  # On the wiki start from Computing/JupyterAtNERSC 
  print('Targets with problematic VI: ', unique_target_csv)
  print('Total number of conflicts to resolve: ', len(unique_targets))

def print_merged_file(vi_gp,output_file):
	vi_gp['Redrock z', 'best z', 'best quality', 'Redrock spectype', 'best spectype', 'all VI issues', 'all VI comments', 'merger comment','N_VI','DELTACHI2', 'ZWARN', 'ZERR','FIBER','FLUX_G', 'FLUX_R', 'FLUX_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MW_TRANSMISSION_G','MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z'].first().to_csv(output_file)

if __name__ == "__main__":
	print('What a cool program.')
