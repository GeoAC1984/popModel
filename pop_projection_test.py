### TO DO: fill value for the resulted projection 

import pandas as pd
import scipy
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from osgeo import gdal
import os, math
from astropy.table import Table
from astropy.convolution import convolve,CustomKernel

path='T:\PopulationModel\Testing\Jamaica'
mask='jampctdl.tif' ##spatial mask; constant
nodata_pop=2147483647
nodata_mask=-3.40282346639e+038
habitable=gdal.Open(os.path.join(path,mask))
mask_arr=habitable.ReadAsArray()
mask_arr = np.array(mask_arr, dtype=np.float64)## assign type 64-bit float to avoiad overload from potentially big numbers
mm = np.ma.masked_values(mask_arr,nodata_mask)# mask no data in spatial mask

scenario='SSP1' ## might be user input; keep as SSP1 for  now
beta=1 
base_year=2020# set the base
proj_table=pd.read_excel(os.path.join(path,'Jamaica_popproj.xlsx'),index_col=[0,1]) ## the table with pop projection numbers
scenario_subset=proj_table.filter(like=scenario,axis=0)## select a subset for scenario

df=pd.read_excel(r'T:\PopulationModel\Testing\Jamaica\distance_window_rd7.xlsx',header=None,sheet='values') ## distance-based weights for Kernel
df=df*13.875 ## roughly converts to 100-m buffer distance  based weights
distance_window=df.as_matrix()
distance_window=np.array(distance_window, dtype=np.float64)
moving_window=np.power(math.e,(-beta*distance_window)) ##multiiplier in the inner parenthesis of the formula
moving_window[moving_window==1]=0 ## reset the 0 for area outside of the buffer

for year in range(len(scenario_subset)):
    print 'base year:', base_year
    projection_year=base_year+10
    in_raster='jampop{}.tif'.format(base_year)

    base_year=projection_year
    np.random.seed(42)
    alpha= np.random.rand(8,18) ## random+ reproducible for now (will be derived from calibration evetially)
    SSP2='jam2010ssp21.tif'
    pop2=scenario_subset['TotalPop_IIASA'].ix[projection_year,scenario]
    print 'projected population number for ', projection_year, 'is ', pop2

    pop=gdal.Open(os.path.join(path,in_raster))
    pop_arr= pop.ReadAsArray()    
    pop_arr = np.array(pop_arr, dtype=np.float64)
    mp = np.ma.masked_values(pop_arr,nodata_pop)
    pop1=np.sum(mp)
    pop_change=pop2-pop1
    
    astropy_kernel=CustomKernel(moving_window)## astropy uses object kernels
    ### ---This is where the potential for each cell is calculated---------------------------------------------------------------------------------------------###
    convolved=np.add(convolve(mp,astropy_kernel),np.power(mp,alpha)) ## inner parenthesis in surface formula + itself in the power of cell-specific alpha
    potential=mm*convolved
    ###--------------------------------------------------------------------------------------------------------------------------------------------------------###    
    if pop_change>0:
        total_potential=potential.sum()
        pot_surface=potential/total_potential
        projected_surface=mp+(pot_surface*abs(pop_change)) ## redistribute pop_change and add to the base year pop
    else:
        potential=1/potential ## reverse rule for likelihood to loose the population
        ## Note:may use different alpha and beta parameters in the future for the years of the population loss; for now will deal with it as above
        total_potential=potential.sum()
        pot_surface=potential/total_potential
        projected_surface=mp-(pot_surface*abs(pop_change)) ## redistribute pop_change and substract from the the base year pop
        
        while (projected_surface<0).any() ==True:## until there are negatvie values in the resulted population distribution, turn negatives to 0 and redistribute the residual
            extra_pop=abs(np.sum(projected_surface[projected_surface<0])) ## get the  pop that needs to be redistributed
            binary_condition=np.where(projected_surface>0,1,0) ## creates a bunary array with 0 in places  where the resulted projecttion=0
            capped_proj_surface=projected_surface.clip(min=0) # capp the projection at 0
            
            '''mask cells in the potential surface with locations of 0 in projected surface
            these cells shouldn't participate in the next iteration of redistributing residual'''
            
            adjusted_pot_surface= binary_condition*pot_surface 
            adjusted_pot_surface=adjusted_pot_surface/adjusted_pot_surface.sum() ## make sure it sums up to 1 
            redistributed_residual=adjusted_pot_surface*extra_pop        
            projected_surface=capped_proj_surface-redistributed_residual#distribute extra pop to the masked surface            
            print projected_surface.sum()
            print 'pop2 difference', projected_surface.sum()-pop2
            
    print 'pop change between base and projection year=', pop_change
    out_projection='jampop{}.tif'.format(projection_year)
    [rows, cols] = pop_arr.shape 
    outfile=os.path.join(path,out_projection)
    driver= gdal.GetDriverByName("GTiff")
    out_raster = driver.Create(outfile, cols, rows, 1, gdal.GDT_Float64)
    proj=pop.GetProjection()
    trans=pop.GetGeoTransform()
    out_raster.GetRasterBand(1).WriteArray(projected_surface)
    out_raster.GetRasterBand(1).SetNoDataValue(np.float64(nodata_pop))
    out_raster.SetGeoTransform(trans)
    out_raster.SetProjection(proj)
    print 'Created projection for ', projection_year

    

