#!/Users/arri/miniconda3/bin/python

import argparse
import os
import time
import numpy as np

# pip install astroquery
from astroquery.mast import Observations

parser = argparse.ArgumentParser(description='2025-05-15 Arno Riffeser (USM@LMU)\npython ./get_TIC_LC.py')
parser.add_argument('numbers', nargs='*', help='TIC number')
parser.add_argument('-prov',  dest='prov', type=str,  default='', help='[%(default)s] provenience name: TESS-SPOC or SPOC or QLP')
parser.add_argument('-name',  dest='name', type=str,  default='', help='[%(default)s] name')
parser.add_argument('-links',  dest='links', action='store_true',  default='', help='[%(default)s] links')
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

if args.prov in ['tess-spoc','spoc','qlp'] :
  provlist=[args.prov]
elif args.prov=='' :
  #provlist=['SPOC','QLP','TESS-SPOC']
  provlist=['tess-spoc','spoc','qlp',]
else :
  print('wrong prov :',args.prov)
  
for n in args.numbers :

  name = 'TIC '+n 
  print(name)
  prename = 'TIC'+n

  for provname in provlist :
  
    print(provname)
    
    obs = Observations.query_criteria(objectname=name, project='TESS',
                                      dataproduct_type='timeseries', provenance_name=provname,
                                      radius='2.5e-08 deg')
    mask = [ url[-7:]=='lc.fits' for url in obs['dataURL'] ]
    obs = obs[mask]
  
    if len(obs)>0 :

      #for i in obs.colnames :
      #  print(i,'  : ',list(obs[i]))

      if args.links :

        obs['dataURL']
        
        for filepath in sorted(obs['dataURL']) :
          fileend=filepath[-7:]
          if fileend=='lc.fits' :
            if provname=='tess-spoc' :
              filesec=filepath[-18:-16]
              filetic=filepath[-32:-22]
              type='t'
            elif provname=='spoc' :
              filesec=filepath[-34:-32]
              filetic=filepath[-25:-15]
              type='s'              
            elif provname=='qlp' :
              filesec=filepath[-37:-35]
              filetic=filepath[-28:-18]
              type='q'              
            fileprov='_'+provname
            if provname=='spoc' and 'fast' in filepath :
              filesec=filepath[-39:-37]
              filetic=filepath[-30:-20]
              type='f'              
              fileprov='_spoc-fast'
            if args.name=='' :
              newname='TIC'+filetic+'_SEC'+filesec+fileprov+'.fits'
            else :
              #newname=args.name+'_TIC'+filetic+'_SEC'+filesec+fileprov+'.fits'
              newname='TOI_'+args.name+'_TIC'+filetic+'_SEC'+filesec+fileprov+'.fits'
            print(newname,'   ','https://mast.stsci.edu/api/v0.1/Download/file/?uri='+filepath)
            # run html list use:
            # awk '(NR>5){printf("./get_TIC_LC.py -links -name \"TOI %04.0f\" %i\n",$2,$3)}' bla | uniq > get_all_TIC_html.sh
            # print('TOI ',args.name,' &#8194; TIC',filetic,'&#8194; SEC',filesec,' &#8194; '+type+' &#8194; <a href="https://mast.stsci.edu/api/v0.1/Download/file/?uri='+filepath+'" download="'+newname+'">'+newname+'</a><br>') 

      else :
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
            for file in sorted(filenames):
              fileend=file[-7:]
              if fileend=='lc.fits' :
                # if provname=='tess-spoc' :
                #   filesec=file[45:47]
                #   filetic=file[31:41]
                # else :
                #   filesec=file[21:23]
                #   filetic=file[30:40]              
                # fileprov='_'+provname
                #
                filepath = os.path.join(dirpath,file)
                if provname=='tess-spoc' :
                  filesec=filepath[-18:-16]
                  filetic=filepath[-32:-22]
                  type='t'
                elif provname=='spoc' :
                  filesec=filepath[-34:-32]
                  filetic=filepath[-25:-15]
                  type='s'              
                elif provname=='qlp' :
                  filesec=filepath[-37:-35]
                  filetic=filepath[-28:-18]
                  type='q'              
                fileprov='_'+provname
                if provname=='spoc' and 'fast' in filepath :
                  filesec=filepath[-39:-37]
                  filetic=filepath[-30:-20]
                  type='f'              
                  fileprov='_spoc-fast'
                if args.name=='' :
                  newname='TIC'+filetic+'_SEC'+filesec+fileprov+'.fits'
                else :
                  newname=args.name+'_TIC'+filetic+'_SEC'+filesec+fileprov+'.fits'
                print('    ',newname)
                os.rename(filepath,newname)

