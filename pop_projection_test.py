### TO DO: create projected_surface fucntion, potential surface function for both positive and negative population growth
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
nodata_ru=65535
nodata_ur=2147483647

habitable=gdal.Open(os.path.join(path,mask))
mask_arr=habitable.ReadAsArray()
mask_arr = np.array(mask_arr, dtype=np.float64)## assign type 64-bit float to avoiad overload from potentially big numbers
mm = np.ma.masked_values(mask_arr,nodata_mask)# mask no data in spatial mask

scenario='SSP1' ## might be user input; keep as SSP1 for  now
urban_beta=1
rural_beta=1
base_year=2000# set the base
proj_table=pd.read_excel(os.path.join(path,'Jamaica_popproj.xlsx'),index_col=[0,1]) ## the table with pop projection numbers
scenario_subset=proj_table.filter(like=scenario,axis=0)## select a subset for scenario

df=pd.read_excel(r'T:\PopulationModel\Testing\Jamaica\distance_window_rd7.xlsx',header=None,sheet='values') ## distance-based weights for Kernel
df=df*13.875 ## roughly converts to 100-m buffer distance  based weights
distance_window=df.as_matrix()
distance_window=np.array(distance_window, dtype=np.float64)
urban_window=np.power(math.e,(-urban_beta*distance_window)) ##multiiplier in the inner parenthesis of the formula
rural_window=np.power(math.e,(-rural_beta*distance_window)) 
urban_window[urban_window==1]=0 ## reset the 0 for area outside of the buffer
rural_window[rural_window==1]=0 

np.random.seed(42)
urban_alpha= np.random.rand(8,18) 
np.random.seed(15)
rural_alpha= np.random.rand(8,18) ## random+ reproducible for now (will be derived from calibration eventially)

for year in range(len(scenario_subset-9)):
    print 'base year:', base_year
    projection_year=base_year+10
    if base_year==2000:
        in_tot_pop='jamtotal{}.tif'.format('001')
        in_urb_pop='jamurban{}.tif'.format('001')
        in_rur_pop='jamrural{}.tif'.format('001')
        
        tot_pop=gdal.Open(os.path.join(path,in_tot_pop))
        tot_pop_arr= tot_pop.ReadAsArray()
        
        urb_pop=gdal.Open(os.path.join(path,in_urb_pop))
        urb_pop_arr= urb_pop.ReadAsArray()
        
        rur_pop=gdal.Open(os.path.join(path,in_rur_pop))
        rur_pop_arr= rur_pop.ReadAsArray() 
    else:
        tot_pop_arr= total_projected_surface
        urb_pop_arr= urban_projected_surface
        rur_pop_arr= rural_projected_surface
        
    base_year=projection_year
    # pop2=scenario_subset['TotalPop_IIASA'].ix[projection_year,scenario]
    pop2_urban=scenario_subset['UrbanPop'].ix[projection_year,scenario]
    pop2_rural=scenario_subset['RuralPop'].ix[projection_year,scenario]
    print 'projected urban population number for ', projection_year, 'is ', pop2_urban
    print 'projected rural population number for ', projection_year, 'is ', pop2_rural

    tot_pop_arr = np.array(tot_pop_arr, dtype=np.float64)
    mtp = np.ma.masked_values(tot_pop_arr,nodata_pop)

    urb_pop_arr = np.array(urb_pop_arr, dtype=np.float64)
    mup = np.ma.masked_values(urb_pop_arr,nodata_ur)

    rur_pop_arr = np.array(rur_pop_arr, dtype=np.float64)
    mrp = np.ma.masked_values(rur_pop_arr,nodata_ru)

    pop1_urban=np.sum(mup)
    urb_pop_change=pop2_urban-pop1_urban

    pop1_rural=np.sum(mrp)
    rur_pop_change=pop2_rural-pop1_rural

    astropy_kernel1=CustomKernel(urban_window)## astropy uses object kernels
    astropy_kernel2=CustomKernel(rural_window)
    ### ---This is where the potential for each cell is calculated---------------------------------------------------------------------------------------------###
    urban_convolved=np.add(convolve(mtp,astropy_kernel1),np.power(mtp,urban_alpha))
    rural_convolved=np.add(convolve(mtp,astropy_kernel2),np.power(mtp,rural_alpha))## inner parenthesis in surface formula + itself in the power of cell-specific alpha
    urban_potential=mm*urban_convolved
    rural_potential=mm*rural_convolved
    print 'pop change between base and projection year=', rur_pop_change+urb_pop_change
    ###--------------------------------------------------------------------------------------------------------------------------------------------------------###    
    if urb_pop_change>0:
        urban_total_potential=urban_potential.sum()
        urban_pot_surface=urban_potential/urban_total_potential
        urban_projected_surface=mup+(urban_pot_surface*abs(urb_pop_change)) ## redistribute pop_change and add to the base year pop

    if rur_pop_change>0:
        rural_total_potential=rural_potential.sum()
        rural_pot_surface=rural_potential/rural_total_potential
        rural_projected_surface=mrp+(rural_pot_surface*abs(rur_pop_change))

    if urb_pop_change<0:
        urban_potential=1/urban_potential ## reverse rule for likelihood to loose the population
        ## Note:may use different alpha and beta parameters in the future for the years of the population loss; for now will deal with it as above
        urban_total_potential=urban_potential.sum()
        urban_pot_surface=urban_potential/urban_total_potential
        urban_projected_surface=mup-(urban_pot_surface*abs(urb_pop_change))## redistribute pop_change and substract from the the base year pop

        while (urban_projected_surface<0).any() ==True:## until there are negatvie values in the resulted population distribution, turn negatives to 0 and redistribute the residual
            print 'capping at 0 and redistributing residual for urban projected surface'
            extra_pop=abs(np.sum(urban_projected_surface[urban_projected_surface<0])) ## get the  pop that needs to be redistributed
            print 'urban extra pop: ',extra_pop
            binary_condition=np.where(urban_projected_surface>0,1,0) ## creates a bunary array with 0 in places  where the resulted projecttion=0
            capped_proj_surface=urban_projected_surface.clip(min=0) # capp the projection at 0
            
            '''mask cells in the potential surface with locations of 0 in projected surface
            these cells shouldn't participate in the next iteration of redistributing residual'''
            
            adjusted_pot_surface= binary_condition*urban_projected_surface 
            adjusted_pot_surface=adjusted_pot_surface/adjusted_pot_surface.sum() ## make sure it sums up to 1 
            redistributed_residual=adjusted_pot_surface*extra_pop        
            urban_projected_surface=capped_proj_surface-redistributed_residual #distribute extra pop to the masked surface            
            print urban_projected_surface.sum()
            print 'pop2 difference', urban_projected_surface.sum()-pop2_urban

    if rur_pop_change<0:
        rural_potential=1/rural_potential 
        rural_total_potential=urban_potential.sum()
        rural_pot_surface=urban_potential/urban_total_potential
        rural_projected_surface=mrp-(rural_pot_surface*abs(rur_pop_change))## redistribute pop_change and substract from the the base year pop        
        
        while (rural_projected_surface<0).any() ==True:## until there are negatvie values in the resulted population distribution, turn negatives to 0 and redistribute the residual
            print 'capping at 0 and redistributing residual for rural projected surface'
            extra_pop=abs(np.sum(rural_projected_surface[rural_projected_surface<0])) ## get the  pop that needs to be redistributed
            print 'rural extra pop: ',extra_pop
            binary_condition=np.where(rural_projected_surface>0,1,0) ## creates a bunary array with 0 in places  where the resulted projecttion=0
            capped_proj_surface=rural_projected_surface.clip(min=0) # capp the projection at 0
            
            '''mask cells in the potential surface with locations of 0 in projected surface
            these cells shouldn't participate in the next iteration of redistributing residual'''
            
            adjusted_pot_surface= binary_condition*rural_projected_surface 
            adjusted_pot_surface=adjusted_pot_surface/adjusted_pot_surface.sum() ## make sure it sums up to 1 
            redistributed_residual=adjusted_pot_surface*extra_pop        
            rural_projected_surface=capped_proj_surface-redistributed_residual #distribute extra pop to the masked surface            
            print rural_projected_surface.sum()
            print 'pop2 difference', rural_projected_surface.sum()-pop2_rural            

    total_projected_surface=urban_projected_surface+rural_projected_surface
    out_projection='jamtotal{}.tif'.format(projection_year)
    out_urban='jamurban{}.tif'.format(projection_year)
    out_rurall='jamrural{}.tif'.format(projection_year)
    [rows, cols] = tot_pop_arr.shape

    outfile=os.path.join(path,out_projection)
    outfile_urb=os.path.join(path,out_urban)
    outfile_rur=os.path.join(path,out_rurall)

##    driver= gdal.GetDriverByName("GTiff")    
##    proj=tot_pop.GetProjection()
##    trans=tot_pop.GetGeoTransform()
##    out_raster0= driver.Create(outfile, cols, rows, 1, gdal.GDT_Float64)
##    out_raster1= driver.Create(outfile_urb, cols, rows, 1, gdal.GDT_Float64)
##    out_raster2= driver.Create(outfile_rur, cols, rows, 1, gdal.GDT_Float64)
##
##    out_raster0.GetRasterBand(1).WriteArray(total_projected_surface)
##    out_raster0.GetRasterBand(1).SetNoDataValue(np.float64(nodata_pop))
##
##    out_raster1.GetRasterBand(1).WriteArray(urban_projected_surface)
##    out_raster1.GetRasterBand(1).SetNoDataValue(np.float64(nodata_pop))
##
##    out_raster2.GetRasterBand(1).WriteArray(rural_projected_surface)
##    out_raster2.GetRasterBand(1).SetNoDataValue(np.float64(nodata_pop))
##
##    out_raster0.SetGeoTransform(trans)
##    out_raster0.SetProjection(proj)
##    out_raster1.SetGeoTransform(trans)
##    out_raster1.SetProjection(proj)
##    out_raster2.SetGeoTransform(trans)
##    out_raster2.SetProjection(proj)

    print 'Created projection for ', projection_year

print 'Done'
