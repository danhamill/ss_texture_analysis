# -*- coding: utf-8 -*-
"""
Created on Mon Jul 03 07:18:29 2017

@author: dan
"""

from __future__ import division
import numpy as np
from sklearn import mixture
from sklearn import cross_validation
from sklearn.metrics import classification_report, confusion_matrix, cohen_kappa_score
from sklearn import preprocessing
import pandas as pd
import os
from osgeo import gdal,ogr,osr
import platform


def read_raster(raster):
    ds = gdal.Open(raster)
    data = ds.GetRasterBand(1).ReadAsArray()
    gt = ds.GetGeoTransform()
    return data, gt

# =========================================================
def get_subs(shp):
    ds = ogr.Open(shp)
    lyr = ds.GetLayer(0)
    a=[]
    for row in lyr:
        a.append(row.substrate)
    lyr.ResetReading()
    del ds
    return a

# =========================================================
def remove_outliers(X, y, k):
   """
   simple outlier removal based on deviation from mean
   """
   mu, sigma = np.mean(X, axis=0), np.std(X, axis=0, ddof=1)
   index = np.all(np.abs((X - mu) / sigma) < k, axis=1)
   return X[index], y[index]

# =========================================================
def CreateRaster(sed_class,gt,outFile):  
    '''
    Exports data to GTiff Raster
    '''
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(26949)
    sed_class = np.squeeze(sed_class)
    sed_class[np.isnan(sed_class)] = -99
    driver = gdal.GetDriverByName('GTiff')
    rows,cols = np.shape(sed_class)
    ds = driver.Create( outFile, cols, rows, 1, gdal.GDT_Float32)      
    if proj is not None:  
        ds.SetProjection(proj.ExportToWkt()) 
    ds.SetGeoTransform(gt)
    ss_band = ds.GetRasterBand(1)
    ss_band.WriteArray(sed_class)
    ss_band.SetNoDataValue(-99)
    ss_band.FlushCache()
    ss_band.ComputeStatistics(False)
    del ds   
   
###################################################################################################################
#####################Entropy, 2 goussicans, 2 sed classes
###################################################################################################################
if __name__ =='__main__':
    if platform.system() == 'Windows':
        clone_root = r'c:\workspace\ss_texture_analysis'
        out_root = clone_root + os.sep + "sedclass_rasters"
        
        
    #input rasters
    ent_dict = {'R01346': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01346_entropy_resampled.tif',
                'R01765': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01765_entropy_resampled.tif',
                'R01767': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01767_entropy_resampled.tif'}
                

    glcm_distributions = clone_root + os.sep + 'glcm_stats' + os.sep + 'merged_GLCM_distributions.csv'
    data = np.genfromtxt(glcm_distributions, delimiter=',', skip_header=1)
    
    data = data[~np.isnan(data).any(axis=1)]
    # for removing outliers
    factor=3 #smaller the number, the more ruthless the pruning
    
    data = data
    data, sedclass = remove_outliers(data[:,:3], data[:,3], factor)
    
    tmp_df = pd.DataFrame({'Entropy':data[:,0],'sedclass':sedclass})
    tmp_df = tmp_df[tmp_df['sedclass'] != 2]
    
    data[:,2] = 1 - data[:,2]
    
    classes = ['Sand','Boulders']
    
    standardize = 0
    
    # split into 50% training, 50% testing
    if standardize==1: # standardize data
       X_train, X_test, y_train, y_test = cross_validation.train_test_split(preprocessing.scale(tmp_df['Entropy'].values), tmp_df['sedclass'].values, test_size=0.5, random_state=0)
       tmp_df['Entropy'] =preprocessing.scale(tmp_df['Entropy'].values)
    else:
       X_train, X_test, y_train, y_test = cross_validation.train_test_split(tmp_df['Entropy'].values, tmp_df['sedclass'].values, test_size=0.5, random_state=0)
    
    #initialize the GMM with means
    g = mixture.GaussianMixture(n_components=2, max_iter=100, random_state=0, covariance_type='tied')
    g.means_init =  np.array([X_train[y_train == i].mean(axis=0) for i in [1,3]])
    g.means_init =np.expand_dims(g.means_init, axis=1) 
    
    # fit the model
    g.fit(np.expand_dims(X_train, axis=1) )
    
     
    #make sure the means are in order
    order = np.argsort(np.squeeze(g.means_))
    g.means_ = g.means_[order]
    try:
       g.covariances_ = g.covariances_[order]
    except:
       pass
    g.weights_ = g.weights_[order]
    
    
    bic = g.bic(np.expand_dims(X_train, axis=1) )
    # test
    y_test_pred = g.predict(np.expand_dims(X_test, axis=1))
    y_test_pred[y_test_pred==1] = 3
    y_test_pred[y_test_pred==0] = 1
    test_accuracy = np.mean(y_test_pred.ravel() == y_test.ravel()) * 100
    print "======================================="
    print "test scores: Entropy"
    print test_accuracy
    
    print(classification_report(y_test_pred.ravel(), y_test.ravel()))
    print (cohen_kappa_score(y_test_pred.ravel(), y_test.ravel()))
    
    # show normalized confusion matrix
    cm = confusion_matrix(y_test.ravel(), y_test_pred.ravel())
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    print(cm)
       
    # Now lets make some maps...
    for (key,v), in zip(ent_dict.items()):         
        ent_raster =v
    
        
        ent_data, gt = read_raster(ent_raster)
    
        
        vec1 = np.c_[ent_data.flatten()]
        vec1[np.isnan(vec1)] = 0; vec1[np.isinf(vec1)] = 0
             
        ind = np.nonzero(vec1)[0]  
        sed_class = np.zeros(np.shape(ent_data.flatten()))*np.nan
        
        for k in xrange(len(ind)):
            sed_class[ind[k]] = g.predict(vec1[ind[k]].reshape(1, -1))
        sed_class[sed_class == 1] = 3
        sed_class[sed_class == 0] = 1
    
        sed_class = np.reshape(sed_class,np.shape(read_raster(ent_raster)[0]))
        sed_class[read_raster(ent_raster)[0]==-99]=np.nan
        
        outFile = out_root + os.sep + str(key) + '_GMM_2class_raster.tif'
    
        CreateRaster(sed_class,gt,outFile)
        
    for (key,v), in zip(ent_dict.items()):         
        ent_raster =v

        ent_data, gt = read_raster(ent_raster)

        vec1 = np.c_[ent_data.flatten()]
        vec1[np.isnan(vec1)] = 0; vec1[np.isinf(vec1)] = 0
             
        ind = np.nonzero(vec1)[0]  
        sed_class = np.zeros(np.shape(ent_data.flatten()))*np.nan
        sed_class = np.c_[sed_class,sed_class]
        for k in xrange(len(ind)):
             sed_class[ind[k]] = g.predict_proba(vec1[ind[k]].reshape(1, -1))[:,g.predict(vec1[ind[k]].reshape(1, -1))[0]]
        sand = sed_class[:,0]
        boulders = sed_class[:,1]
        
        sand, boulders =[np.reshape(i,np.shape(read_raster(ent_raster)[0]))for i in [sand, boulders]]
        
        sand[read_raster(ent_raster)[0]==-99]=np.nan
        boulders[read_raster(ent_raster)[0]==-99]=np.nan
        
        outFile = out_root + os.sep + str(key) + '_GMM2_sand_proba.tif'
        CreateRaster(sand,gt,outFile)
        
        outFile = out_root + os.sep + str(key) + '_GMM2_boulders_proba.tif'
        CreateRaster(boulders,gt,outFile)