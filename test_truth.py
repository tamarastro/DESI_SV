import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table, join, vstack

tiledir   = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/'
tiles = ['68002']
nights = ['20200315']
petals = ['0','1', '2', '3', '4', '5', '6' ,'7', '8', '9']

directory_truth_tables = '/global/project/projectdirs/desi/sv/vi/TruthTables/'
truth_table = 'truth_table_QSO_v1.2.csv'

vi = pd.read_csv(directory_truth_tables+truth_table, delimiter = ",", engine='python')

tspec = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[0])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='ZBEST')
for i in range(1,len(petals)):
    tspecn = Table.read(tiledir+'/'+tiles[0] + '/'+nights[0]+'/zbest-'+str(petals[i])+'-'+str(tiles[0])+'-'+nights[0]+'.fits',hdu='ZBEST')
    tspec = vstack([tspec,tspecn])
tspec_check = tspec.to_pandas()
print(tspec_check.keys())

tspec_df = tspec['TARGETID','CHI2'].to_pandas()

vi = vi.merge(tspec_df, how='left', on='TARGETID')
print(vi.keys())
