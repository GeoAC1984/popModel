import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import os, math

nodata_pop=2147483647
nodata_mask=-3.40282346639e+038
nodata_ru=65535
nodata_ur=2147483647

def read_tif(path,tif,nodata):
    pop=gdal.Open(os.path.join(path,tif))
    pop_arr= pop.ReadAsArray()
    pop_arr = np.array(pop_arr, dtype=np.float64)
    mp = np.ma.masked_values(pop_arr,nodata)
    return mp

path=r'T:\PopulationModel\Testing\Jamaica'
projections=r'T:\PopulationModel\Testing\Jamaica\projection'
## readinng in observed

tot_pop='jamtotal{}.tif'.format('001')
urb_pop='jamurban{}.tif'.format('001')
rur_pop='jamrural{}.tif'.format('001')

mod_tot_pop='jamtotal{}.tif'.format('2000')
mod_urb_pop='jamurban{}.tif'.format('2000')
mod_rur_pop='jamrural{}.tif'.format('2000')

mtp_obs = read_tif(path,tot_pop,nodata_pop)
mup_obs =read_tif(path,urb_pop,nodata_ur)
mrp_obs = read_tif(path,rur_pop,nodata_ru)

mtp_mod = read_tif(projections,mod_tot_pop,nodata_pop)
mup_mod =read_tif(projections,mod_urb_pop,nodata_pop)
mrp_mod= read_tif(projections,mod_rur_pop,nodata_pop)

