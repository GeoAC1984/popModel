import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import os, math
from scipy.optimize import curve_fit

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

path=r'C:\GravityModel\Jamaica'

## for this part -- use 1990 (actual) as observed and 2000 (actual) as modeled
tot_pop='jamtotal{}.tif'.format('90')
urb_pop='jamurban{}.tif'.format('90')
rur_pop='jamrural{}.tif'.format('90')

mod_tot_pop='jamtotal{}.tif'.format('00')
mod_urb_pop='jamurban{}.tif'.format('00')
mod_rur_pop='jamrural{}.tif'.format('00')

mtp_obs = read_tif(path,tot_pop,nodata_pop)
mup_obs =read_tif(path,urb_pop,nodata_pop)
mrp_obs = read_tif(path,rur_pop,nodata_pop)

mtp_mod = read_tif(path,mod_tot_pop,nodata_pop)
mup_mod =read_tif(path,mod_urb_pop,nodata_ur)
mrp_mod= read_tif(path,mod_rur_pop,nodata_ru)


## find B that will minimize Err
ydata=mtp_mod/mtp_mod.sum()
y=ydata.flatten()
xdata=mtp_obs/mtp_obs.sum()
x=xdata.flatten()
### assume obs (2000)  as first-order polynomial function, where c is 1990, x are errors and m is beta param
### obs=c+mx
### p_coeff will return coefficients for c (1) and beta param
p_coeff= np.polynomial.polynomial.polyfit(x, y, 1)
print(p_coeff) 
