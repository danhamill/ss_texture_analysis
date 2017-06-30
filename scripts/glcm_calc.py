# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:12:39 2017

@author: dan
"""

import gdal, osr
import numpy as np
from scipy.interpolate import RectBivariateSpline
from numpy.lib.stride_tricks import as_strided as ast
import dask.array as da
from joblib import Parallel, delayed
import os
from skimage.feature import greycomatrix, greycoprops
import subprocess
import platform

def im_resize(im,Nx,Ny):
    '''
    resize array by bivariate spline interpolation
    '''
    ny, nx = np.shape(im)
    xx = np.linspace(0,nx,Nx)
    yy = np.linspace(0,ny,Ny)
    
    try:
        im = da.from_array(im, chunks=1000)   #dask implementation
    except:
        pass

    newKernel = RectBivariateSpline(np.r_[:ny],np.r_[:nx],im)
    return newKernel(yy,xx)
    
def mean_var(P):
    '''
    Function go calculate GLCM mean and GLCM Var directly from a GLCM
    '''
    (num_level, num_level2, num_dist, num_angle) = P.shape
    I, J = np.ogrid[0:num_level, 0:num_level]
    I = np.array(range(num_level)).reshape((num_level, 1, 1, 1))
    mean_i = np.apply_over_axes(np.sum, (I * P), axes=(0, 1))[0, 0]
    diff_i = I - np.apply_over_axes(np.sum, (I * P), axes=(0, 1))[0, 0]
    var_h = np.apply_over_axes(np.sum, (P * (diff_i) ** 2), axes=(0, 1))[0, 0]
    return mean_i, var_h 
    
def entropy_calc(glcm):
    '''
    Funciton go calculate entropy from a GLCM
    '''
    with np.errstate(divide='ignore', invalid='ignore'):
        horizontal_entropy = np.apply_over_axes(np.nansum,(np.log(glcm)*-1*glcm),axes=(0,1))[0,0]
        horizontal_entropy = np.asarray([[horizontal_entropy[0,0]]])
    return horizontal_entropy
    
def p_me(Z, win,dist,angle):
    '''
    Function to calculate GLCM and subsuquent GLCM properties
    '''
    if np.count_nonzero(Z) > 0.75*win**2: 
        glcm = greycomatrix(Z, [dist], [angle], 256, symmetric=True, normed=True)
        homo = greycoprops(glcm, 'homogeneity')
        entropy = entropy_calc(glcm)
        mean, var = mean_var(glcm)
        return (homo, entropy, var)
    else:
        return (0,0,0)
        
        
        
def read_raster(in_raster):
    '''
    Function to read GTiff rasters
    Returns: data, Easting grid, northing grid, and geotransform
    '''
    ds = gdal.Open(in_raster)
    data = ds.GetRasterBand(1).ReadAsArray()
    data[data<=0] = np.nan
    gt = ds.GetGeoTransform()
    xres = gt[1]
    yres = gt[5]
    
    # get the edge coordinates and add half the resolution 
    # to go to center coordinates
    xmin = gt[0] + xres * 0.5
    xmax = gt[0] + (xres * ds.RasterXSize) - xres * 0.5
    ymin = gt[3] + (yres * ds.RasterYSize) + yres * 0.5
    ymax = gt[3] - yres * 0.5
    del ds
    # create a grid of xy coordinates in the original projection
    xx, yy = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
    return data, xx, yy, gt
    
def norm_shape(shap):
   '''
   Normalize numpy array shapes so they're always expressed as a tuple,
   even for one-dimensional shapes.
   '''
   try:
      i = int(shap)
      return (i,)
   except TypeError:
      # shape was not a number
      pass

   try:
      t = tuple(shap)
      return t
   except TypeError:
      # shape was not iterable
      pass

   raise TypeError('shape must be an int, or a tuple of ints')

def sliding_window(a, ws, ss = None, flatten = True):
    '''
    Return a sliding window over a in any number of dimensions
    '''      
    if None is ss:
        # ss was not provided. the windows will not overlap in any direction.
        ss = ws
    ws = norm_shape(ws)
    ss = norm_shape(ss)
    # convert ws, ss, and a.shape to numpy arrays
    ws = np.array(ws)
    ss = np.array(ss)
    shap = np.array(a.shape)
    # ensure that ws, ss, and a.shape all have the same number of dimensions
    ls = [len(shap),len(ws),len(ss)]
    if 1 != len(set(ls)):
        raise ValueError(\
        'a.shape, ws and ss must all have the same length. They were %s' % str(ls))
    
    # ensure that ws is smaller than a in every dimension
    if np.any(ws > shap):
        raise ValueError(\
        'ws cannot be larger than a in any dimension.\
     a.shape was %s and ws was %s' % (str(a.shape),str(ws)))
     
    # how many slices will there be in each dimension?
    newshape = norm_shape(((shap - ws) // ss) + 1)
    
    
    # the shape of the strided array will be the number of slices in each dimension
    # plus the shape of the window (tuple addition)
    newshape += norm_shape(ws)
    
    
    # the strides tuple will be the array's strides multiplied by step size, plus
    # the array's strides (tuple addition)
    newstrides = norm_shape(np.array(a.strides) * ss) + a.strides
    a = ast(a,shape = newshape,strides = newstrides)
    if not flatten:
        return a
    # Collapse strided so that it has one more dimension than the window.  I.e.,
    # the new array is a flat list of slices.
    meat = len(ws) if ws.shape else 0
    firstdim = (np.product(newshape[:-meat]),) if ws.shape else ()
    dim = firstdim + (newshape[-meat:])
    # remove any dimensions with size 1
    dim = filter(lambda i : i != 1,dim)
    
    return a.reshape(dim), newshape
    
def CreateRaster(data,gt,proj,outFile):  
    '''
    Function to write numpy array to GTiff rasters
    '''
    data = np.squeeze(data)
    data[np.isnan(data)] = -99
    data[data>100] = -99
    driver = gdal.GetDriverByName('GTiff')
    rows,cols = np.shape(data)
    ds = driver.Create( outFile, cols, rows, 1, gdal.GDT_Float32)      
    if proj is not None:  
        ds.SetProjection(proj.ExportToWkt()) 
    ds.SetGeoTransform(gt)
    band = ds.GetRasterBand(1)
    band.WriteArray(data)
    band.SetNoDataValue(-99)
    band.FlushCache()
    band.ComputeStatistics(False)
    del ds
    
    
if __name__ == '__main__':  
    
    #Set up root directories for input files and output files
    if platform.system() == 'Windows':
        clone_root = r'c:\workspace\ss_texture_analysis'
        out_root = clone_root + os.sep + "glcm_rasters"
    
    epsg_code=26949
    win = 12
    
    #Build paths to input side scan sonar rasters
    ss_dict = {'R01346': clone_root + os.sep + 'ss_rasters' + os.sep + 'R01346.tif',
               'R01765': clone_root + os.sep + 'ss_rasters' + os.sep + 'R01765.tif',
               'R01767': clone_root + os.sep + 'ss_rasters' + os.sep + 'R01767.tif'}

    
    for (k,v) in ss_dict.items():
        hom_file = out_root + os.sep + k + "_homo.tif"
        ent_file = out_root + os.sep + k +  "_entropy.tif"
        var_file = out_root + os.sep + k +  "_var.tif"

        #Read side scan sonar raster and prep for GLCM calulations
        merge, xx, yy, gt = read_raster(v)
        Ny, Nx = np.shape(merge)
        merge[np.isnan(merge)] = 0
        
        #Calculate sliding window with no overlap in either direction
        Z,ind = sliding_window(merge,(win,win),(win,win))
        
        #Calculte GLCM
        w = Parallel(n_jobs = 1, verbose=0)(delayed(p_me)(Z[k], win, 5, 0) for k in xrange(len(Z)))
        homo, ent, var = zip(*w)
    
        
        # Reshape GLCM texture features into array with sane number of windows
        plt_homo = np.reshape(homo , ( ind[0], ind[1] ) )
        plt_ent = np.reshape(ent, (ind[0],ind[1]))
        plt_var = np.reshape(var, (ind[0],ind[1]))    
        del homo, ent, var
        
        #Resize Images to input raster cell size
        homogeneity = im_resize(plt_homo,Nx,Ny)
        homogeneity[merge==0]=np.nan
        homogeneity[homogeneity<0]=np.nan
    
        ENT = im_resize(plt_ent,Nx,Ny)
        ENT[merge==0]=np.nan
        ENT[ENT<0]=np.nan
    
        var = im_resize(plt_var,Nx,Ny)
        var[merge==0]=np.nan
        var[var<0] = np.nan
        del plt_homo, plt_ent, plt_var, w, Z, ind, Ny, Nx
    
        
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(epsg_code)
        
        #Write GLCM porperty arrays to GTiff raster
        CreateRaster( homogeneity, gt, proj,hom_file)
        CreateRaster( ENT, gt, proj,ent_file)
        CreateRaster(var, gt, proj,var_file)
        
        #Resampled file names
        hom_resamp = hom_file[0:-4] + '_resampled.tif'
        ent_resamp = ent_file[0:-4] + '_resampled.tif'
        var_resamp = var_file[0:-4] + '_resampled.tif'
        
        #Resample GLCM properties to the window size
        subprocess.call(['gdalwarp','-tr', str(win/4), str(win/4), '-r', 'average', hom_file, hom_resamp] )
        subprocess.call(['gdalwarp','-tr', str(win/4), str(win/4), '-r', 'average', ent_file, ent_resamp] )
        subprocess.call(['gdalwarp','-tr', str(win/4), str(win/4), '-r', 'average', var_file, var_resamp] )
        
        del ent_file, hom_file, var_file, ent_resamp, hom_resamp, var_resamp, merge, k, v, xx, yy, ENT, homogeneity, var, gt, proj
    