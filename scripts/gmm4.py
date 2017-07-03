# -*- coding: utf-8 -*-
"""
Created on Mon Jul 03 08:17:33 2017

@author: dan
"""

from __future__ import division
import numpy as np
from sklearn import mixture
from sklearn.metrics import classification_report, confusion_matrix
import pandas as pd
import os
from osgeo import gdal,ogr,osr
import platform
from sklearn.model_selection import train_test_split

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
    hom_dict = {'R01346': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01346_homo_resampled.tif',
                'R01765': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01765_homo_resampled.tif',
                'R01767': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01767_homo_resampled.tif'}  
    
    var_dict = {'R01346': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01346_var_resampled.tif',
                'R01765': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01765_var_resampled.tif',
                'R01767': clone_root + os.sep + 'glcm_rasters' + os.sep + 'R01767_var_resampled.tif'}  

    classes = ['Sand','Gravel','Boulders']

    glcm_distributions = clone_root + os.sep + 'glcm_stats' + os.sep + 'merged_GLCM_distributions.csv'
    data = np.genfromtxt(glcm_distributions, delimiter=',', skip_header=1)
    
    data = data[~np.isnan(data).any(axis=1)]
    # for removing outliers
    factor=3 #smaller the number, the more ruthless the pruning
    
    
    data = data
    data, sedclass = remove_outliers(data[:,:3], data[:,3], factor)
    data[:,2] = 1 - data[:,2]
    tmp_df = pd.DataFrame({'Variance':data[:,1],'Homogeneity':data[:,2],'sedclass':sedclass})



    
    X_train, X_test, y_train, y_test = train_test_split(tmp_df[['Variance','Homogeneity']].values, tmp_df['sedclass'].values, test_size=0.5, random_state=0)

         
    #initialize the GMM with means
    g = mixture.GaussianMixture(n_components=4, max_iter=100, random_state=0, covariance_type='full')
      
    means =  np.array([X_train[y_train == i].mean(axis=0) for i in range(1,len(classes)+1)])
      
    df = pd.DataFrame(data=means,columns=['Variance','Homogeneity'])
      
    df2 = pd.DataFrame(index=['sand','sand_gravel','gravel_boluders','boulders'],columns = ['Variance','Homogeneity'])
    df2.loc['sand'] = pd.Series({'Variance':df.iloc[int(0)].Variance,
                                    'Homogeneity':df.iloc[int(0)].Homogeneity})
    df2.loc['sand_gravel'] = pd.Series({'Variance': (df.iloc[1].Variance+ df.iloc[0].Variance)/2,'Homogeneity': (df.iloc[1].Homogeneity+ df.iloc[0].Homogeneity)/2})
    df2.loc['gravel_boluders'] = pd.Series({'Variance': (df.iloc[2].Variance+ df.iloc[1].Variance)/2,'Homogeneity': 
                                          (df.iloc[2].Homogeneity+ df.iloc[1].Homogeneity)/2})
    df2.loc['boulders'] = pd.Series({'Variance':df.iloc[2].Variance,
                                      'Homogeneity':df.iloc[2].Homogeneity})
    means = df2.values
    del df, df2
    means = means.astype('float')
      
    g.means_init = means
    
    # fit the model
    g.fit(X_train )
   

        
    bic = g.bic(X_train)

    # test
    y_test_pred = g.predict(X_test)
    y_test_pred[y_test_pred==1]=2
    y_test_pred[y_test_pred==0]=1
    test_accuracy = np.mean(y_test_pred.ravel() == y_test.ravel()) * 100
    print "======================================="
    print "test scores: Homogeneity and GLCM Variance"
    print test_accuracy

    print(classification_report(y_test_pred.ravel(), y_test.ravel()))

    print "======================================="
    print "Normalized Confusion Matrix:"  
    cm = confusion_matrix(y_test.ravel(), y_test_pred.ravel())
    cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    print(cm)

    for (key,v), (key1,v1)in zip(var_dict.items(), hom_dict.items())[0:1]:         
        var_raster =v
        homo_raster =v1
        
        var_data, gt = read_raster(var_raster)
        homo_data = read_raster(homo_raster)[0]
        
        vec1 = np.c_[var_data.flatten(),1-homo_data.flatten()]
        vec1[np.isnan(vec1)] = 0; vec1[np.isinf(vec1)] = 0
             
        ind = np.nonzero(vec1)[0]  
        sed_class = np.zeros(np.shape(var_data.flatten()))*np.nan
        
        for k in xrange(len(ind)):
            sed_class[ind[k]] = g.predict(vec1[ind[k]].reshape(1, -1))
        sed_class[sed_class == 1] = 2
        sed_class[sed_class == 0] = 1
    
        sed_class = np.reshape(sed_class,np.shape(read_raster(homo_raster)[0]))
        sed_class[read_raster(var_raster)[0]==-99]=np.nan
        
        outFile = out_root + os.sep + str(key) + '_GMM4_sedclass.tif'
        CreateRaster(sed_class,gt,outFile)
    
    for (key,v), (key1,v1)in zip(var_dict.items(), hom_dict.items())[0:1]:         
        var_raster =v
        homo_raster =v1
        
        var_data, gt = read_raster(var_raster)
        homo_data = read_raster(homo_raster)[0]
        
        vec1 = np.c_[var_data.flatten(),1-homo_data.flatten()]
        vec1[np.isnan(vec1)] = 0; vec1[np.isinf(vec1)] = 0
             
        ind = np.nonzero(vec1)[0]  
        sed_class = np.zeros(np.shape(var_data.flatten()))*np.nan
        
        sed_class = np.c_[sed_class,sed_class,sed_class,sed_class]
        
        for k in xrange(len(ind)):
             sed_class[ind[k]] = g.predict_proba(vec1[ind[k]].reshape(1, -1))
         
        sand = sed_class[:,0] 
        gravel = np.add(sed_class[:,1],sed_class[:,2])
        boulders = sed_class[:,3] 
        sand, gravel, boulders =[np.reshape(i,np.shape(read_raster(homo_raster)[0]))for i in [sand, gravel, boulders]]
        sand[read_raster(var_raster)[0]==-99]=np.nan
        gravel[read_raster(var_raster)[0]==-99]=np.nan
        boulders[read_raster(var_raster)[0]==-99]=np.nan
        
        outFile = out_root + os.sep + str(key) + '_GMM4_sand_proba.tif'
        CreateRaster(sand,gt,outFile)
        outFile = out_root + os.sep + str(key) + '_GMM4_gravel_proba.tif'
        CreateRaster(gravel,gt,outFile)          
        outFile = out_root + os.sep + str(key) + '_GMM4_boulders_proba.tif'
        CreateRaster(boulders,gt,outFile)