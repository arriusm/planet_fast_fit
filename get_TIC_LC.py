#!/Users/arri/miniconda3/bin/python

import argparse
import os
import time
import numpy as np

# pip install astroquery
from astroquery.mast import Observations

parser = argparse.ArgumentParser(description='2025-05-09 Arno Riffeser (USM@LMU)\npython ./get_TIC_LC.py')
parser.add_argument('numbers', nargs='*', help='TIC number')
parser.add_argument('-prov',  dest='prov', type=str,  default='SPOC', help='[%(default)s] provenience name: TESS-SPOC or SPOC or QLP')
args = parser.parse_args()

# https://archive.stsci.edu/hlsp/qlp/provenancesearch
# 
# In Portal
# 
# Mission = "TESS" for mission-produced light curves, "HLSP" for High Level Science Products.
#           If you want to only include official, mission-produced data, and exclude all HLSP products, set MISSION = "TESS".
# 
# Project = "TESS" can be used to find any TESS-related products.
#           This includes official, mission-produced data AND High Level Science Products
#           that are derived from or relate to TESS data.
#           If you want the broadest possible search for "anything related to TESS", set PROJECT = "TESS".
# 
# Provenance = "QLP" can be used to find any QLP light curves.
#              This will exclude all official, mission-produced light curves,
#              and any other HLSP that are not from the "QLP" HLSP.
#              Each HLSP has the Provenance name set to their HLSP shortname, and of
#              course you can include more than one in a given search.
#
# see
#  import lightkurve as lk
#  lk.search_lightcurve("TIC 40292751",author = 'TESS-SPOC')[0].download(flux_column="pdcsap_flux")


for n in args.numbers :

  name = 'TIC '+n 
  prename = 'TIC'+n  
  provname=args.prov

  # e.g.
  # name = 'TIC 407905370' # TOI 3573
  # name = 'TIC 468692991' # TOI 6136
  # name = 'TIC 155258153' # TOI 6195
  # name = 'TIC 15682927'  # TOI 5886
  # provname='QLP'
  
  print(name)
  print(provname)

  obs = Observations.query_criteria(objectname=name, project='TESS', dataproduct_type='timeseries', provenance_name=provname, radius='2.5e-08 deg')
  mask = [ url[-7:]=='lc.fits' for url in obs['dataURL'] ]
  obs = obs[mask]
  
  for i in obs['dataURL'] :
      print('https://mast.stsci.edu/api/v0.1/Download/file/?uri='+i)

  #for i in obs.colnames :
  #  print(i,'  : ',list(obs[i]))
    
  if len(obs)>0 :
  
    data_products = Observations.get_product_list(obs)
    mask = [ url[-7:]=='lc.fits' for url in data_products['productFilename'] ]
    data_products = data_products[mask]
  
    #for i in data_products.colnames :
    #  print(i,'  : ',list(data_products[i]))

    manifest = Observations.download_products(data_products, productType="SCIENCE")

    timestamp = time.strftime('%y%m%d%H%M%S')
    downdir =  'mastDownload_'+timestamp
    os.rename('mastDownload',downdir)

    walker = os.walk(downdir)
    for dirpath, dirnames, filenames in walker:
      if filenames:
        for file in filenames:
          fileend=file[-7:]
          if fileend=='lc.fits' :
            filesec=file[21:23]
            filetic=file[30:40]
            filepath = os.path.join(dirpath,file)
            newname='TIC'+filetic+'_SEC'+filesec+'_'+provname+'.fits'
            print('    ',newname)
            os.rename(filepath,newname)

