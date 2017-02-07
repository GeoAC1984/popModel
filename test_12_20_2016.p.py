import pandas as pd
import scipy
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import os, math
from astropy.table import Table
from astropy.convolution import convolve,CustomKernel

base_year=2050 # set the base
projection_year=base_year+10
beta=1 
alpha=1
nodata_pop=2147483647
nodata_mask=-3.40282346639e+038

path='T:\PopulationModel\Testing\Jamaica'
raster='jampop.tif'
mask='jampctdl.tif'
SSP2='jam2010ssp21.tif'

proj_table=pd.read_excel(os.path.join(path,'Jamaica_popproj.xlsx'),index_col=[0,1])
pop_slice= proj_table.xs([projection_year,'SSP1'],level=['Year','Scenario'])
pop2=pop_slice.iloc[0]['TotalPop_IIASA']

df=pd.read_excel(r'C:\GravityModel\popproj_oct2016\distance_window_rd7.xlsx',header=None,sheet='values')
df=df*13.875
window=df.as_matrix()
window=np.array(window, dtype=np.float64)
moving_window=np.power(math.e,(-beta*window))
moving_window[moving_window==1]=0 ## anything in the power of zero gets set to 1==>set it to zero

pop=gdal.Open(os.path.join(path,raster))
habitable=gdal.Open(os.path.join(path,mask))
ssp2_pop=gdal.Open(os.path.join(path,SSP2))

pop_arr= pop.ReadAsArray()
mask_arr=habitable.ReadAsArray()
ssp2_arr=ssp2_pop.ReadAsArray()

mask_arr = np.array(mask_arr, dtype=np.float64)## assign type 64-bit float to avoiad overload from potentially big numbers
pop_arr = np.array(pop_arr, dtype=np.float64)
ssp2_arr = np.array(ssp2_arr, dtype=np.float64)

mp = np.ma.masked_values(pop_arr,nodata_pop)
base_pop=np.sum(mp)
pop_change=pop2-base_pop


##pop_arr[pop_arr==nodata_pop]=0## turn nodata into 0; only needed for scipy convolution
[rows, cols] = pop_arr.shape 

mm = np.ma.masked_values(mask_arr,nodata_mask)# mask no data in spatial mask
ssp_m=np.ma.masked_values(ssp2_arr,nodata_pop)# mask no data in ssp2 model

## Decided to go with astropy's convolve for ability to handle NoData
# potential=ndimage.filters.convolve(pop_arr,moving_window,mode='nearest')## scipy convolution
astropy_kernel=CustomKernel(moving_window)## astropy uses object kernels
# potential=ndimage.filters.convolve(pop_arr,moving_window,mode='nearest')+pop_arr))
convolved=np.add(convolve(mp,astropy_kernel),mp)## astropy convolution; inner parenthesis in surface formula + itself
if pop_change>0:
    potential=mm*convolved ## outer part of surface the formula
else:
    potential=1(mm*convolved)
total_potential=potential.sum()
pot_surface=potential/total_potential## redistribute into surface
projected_surface=pot_surface*pop_change+mp ## redistribute pop_change and add to the base year pop


### TO DO:
## if pop_change<0
## inverse the surface
## cap the projected surface at 0


# some diagnostics
errors=ssp_m-projected_surface
mean_e=np.mean(errors)
max_e=np.max(errors)
min_e=np.min(errors)

## fill no data with 0 to write out values in raster
errors_tif=np.ma.filled(errors,fill_value=0)
errors_tif2=np.asarray(errors_tif,dtype=np.int64)
pop_downscaled=np.ma.filled(projected_surface,fill_value=0)
pop_downscaled=np.asarray(pop_downscaled,dtype=np.int64)
##
##
####outfile=r'T:\PopulationModel\Testing\Jamaica\projected_surface.tif'
####driver= gdal.GetDriverByName("GTiff")
####out_raster = driver.Create(outfile, cols, rows, 1, gdal.GDT_Float64)
####proj=pop.GetProjection()
####trans=pop.GetGeoTransform()
####out_raster.GetRasterBand(1).WriteArray(pop3)
####out_raster.GetRasterBand(1).SetNoDataValue(0)
####out_raster.SetGeoTransform(trans)
####out_raster.SetProjection(proj)
####
####outfile2=r'T:\PopulationModel\Testing\Jamaica\errors.tif'
####driver= gdal.GetDriverByName("GTiff")
####out_raster2 = driver.Create(outfile2, cols, rows, 1, gdal.GDT_Float64)
####proj2=pop.GetProjection()
####trans2=pop.GetGeoTransform()
####out_raster2.GetRasterBand(1).WriteArray(errors_tif2)
####out_raster2.GetRasterBand(1).SetNoDataValue(0)
####out_raster2.SetGeoTransform(trans2)
####out_raster2.SetProjection(proj2)
