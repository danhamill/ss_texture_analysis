# -*- coding: utf-8 -*-
"""
Created on Sun Jul 02 14:48:22 2017

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
csv = r"C:\workspace\GLCM\new_output\merged_aggregraded_distributions.csv"
data = np.genfromtxt(csv, delimiter=',', skip_header=1)

data = data[~np.isnan(data).any(axis=1)]
# for removing outliers
factor=3 #smaller the number, the more ruthless the pruning

data = data
data, sedclass = remove_outliers(data[:,:3], data[:,3], factor)

tmp_df = pd.DataFrame({'Entropy':data[:,0],'sedclass':sedclass})
tmp_df = tmp_df[tmp_df['sedclass'] != 2]

data[:,2] = 1 - data[:,2]
predictors = ['Entropy','Variance','Homogeneity']



classes = ['Sand','Boulders']

standardize = 0

for covtype in ['tied']:
   print 'Working on covariance type %s...' %(covtype,)
   
   for n in [0]:
      print "working on "+predictors[n]
      print "Working on GMM..."
      # split into 50% training, 50% testing
      if standardize==1: # standardize data
         X_train, X_test, y_train, y_test = cross_validation.train_test_split(preprocessing.scale(tmp_df['Entropy'].values), tmp_df['sedclass'].values, test_size=0.5, random_state=0)
         tmp_df['Entropy'] =preprocessing.scale(tmp_df['Entropy'].values)
      else:
         X_train, X_test, y_train, y_test = cross_validation.train_test_split(tmp_df['Entropy'].values, tmp_df['sedclass'].values, test_size=0.5, random_state=0)
    
      #initialize the GMM with means
      g = mixture.GaussianMixture(n_components=2, max_iter=100, random_state=0, covariance_type=covtype)
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
      print "test scores: "+predictors[n]
      print test_accuracy

      print(classification_report(y_test_pred.ravel(), y_test.ravel()))
      print (cohen_kappa_score(y_test_pred.ravel(), y_test.ravel()))

      # show normalized confusion matrix
      cm = confusion_matrix(y_test.ravel(), y_test_pred.ravel())
      cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
      print(cm)
         
      # Now lets make some maps...
      ent_dict = {'R01346':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_04\3\R01346_R01347_3_entropy_resampled.tif",
                'R01765':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_09_2\3\R01765_3_entropy_resampled.tif",
                'R01767':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_09\3\R01767_3_entropy_resampled.tif"}      
      rasters = []
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
          
          outFile = r'C:\workspace\GLCM\new_output\gmm_rasters' + os.sep + str(key) + '_GMM_2class_raster.tif'
          rasters.append(outFile)
          CreateRaster(sed_class,gt,outFile)
          
      for (key,v), in zip(ent_dict.items()):         
          ent_raster =v

          
          ent_data, gt = read_raster(ent_raster)

          
          vec1 = np.c_[ent_data.flatten()]
          vec1[np.isnan(vec1)] = 0; vec1[np.isinf(vec1)] = 0
               
          ind = np.nonzero(vec1)[0]  
          sed_class = np.zeros(np.shape(ent_data.flatten()))*np.nan
          
          for k in xrange(len(ind)):
               sed_class[ind[k]] = g.predict_proba(vec1[ind[k]].reshape(1, -1))[:,g.predict(vec1[ind[k]].reshape(1, -1))[0]]
               
#          sed_class[sed_class == 1] = 2
#          sed_class[sed_class == 0] = 1
      
          sed_class = np.reshape(sed_class,np.shape(read_raster(ent_raster)[0]))
          
          outFile = r'C:\workspace\GLCM\new_output\gmm_rasters' + os.sep + str(key) + '_GMM2class_proba_raster.tif'
          CreateRaster(sed_class,gt,outFile)


      
###################################################################################################################
#####################Variance + Homogeneity,4 goussicans
###################################################################################################################
csv = r"C:\workspace\GLCM\new_output\merged_aggregraded_distributions.csv"
data = np.genfromtxt(csv, delimiter=',', skip_header=1)

data = data[~np.isnan(data).any(axis=1)]
# for removing outliers
factor=3 #smaller the number, the more ruthless the pruning

data = data
data, sedclass = remove_outliers(data[:,:3], data[:,3], factor)
data[:,2] = 1 - data[:,2]
tmp_df = pd.DataFrame({'Variance':data[:,1],'Homogeneity':data[:,2],'sedclass':sedclass})



predictors = ['Entropy','Variance','Homogeneity']



classes = ['Sand','Gravel','Boulders']

standardize = 0

for covtype in ['full']:
   print 'Working on covariance type %s...' %(covtype,)
   
   for n in [1]:
      print "working on "+predictors[n]
      print "Working on GMM..."
      # split into 50% training, 50% testing
      if standardize==1: # standardize data
         X_train, X_test, y_train, y_test = cross_validation.train_test_split(preprocessing.robust_scale(tmp_df[['Variance','Homogeneity']].values), tmp_df['sedclass'].values, test_size=0.5, random_state=0)
         tmp_df['Variance'] = preprocessing.robust_scale(tmp_df[['Variance']].values)
         tmp_df['Homogeneity'] = preprocessing.robust_scale(tmp_df[['Variance']].values)        
      else:
         X_train, X_test, y_train, y_test = cross_validation.train_test_split(tmp_df[['Variance','Homogeneity']].values, tmp_df['sedclass'].values, test_size=0.5, random_state=0)

         
      #initialize the GMM with means
      g = mixture.GaussianMixture(n_components=4, max_iter=100, random_state=0, covariance_type=covtype)
      
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
      
      #g.means_init =  np.array([X_train[y_train == i].mean(axis=0) for i in range(1,len(classes)+1)])
      g.means_init = means
      # fit the model
      g.fit(X_train )
   
      #make sure the means are in order
      order = np.argsort(g.means_[:,0])
      g.means_[:,0] = g.means_[:,0][order]
      g.covariances_[:,0] = g.covariances_[:,0][order]
        
      bic = g.bic(X_train)

      # test
      y_test_pred = g.predict(X_test)
      y_test_pred[y_test_pred==1]=2
      y_test_pred[y_test_pred==0]=1
      test_accuracy = np.mean(y_test_pred.ravel() == y_test.ravel()) * 100
      print "======================================="
      print "test scores: "+predictors[n]
      print test_accuracy

      print(classification_report(y_test_pred.ravel(), y_test.ravel()))
      print (cohen_kappa_score(y_test_pred.ravel(), y_test.ravel()))
      # show normalized confusion matrix
      cm = confusion_matrix(y_test.ravel(), y_test_pred.ravel())
      cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
      print(cm)
        
      
      # Now lets make some maps...
      var_dict = {'R01346':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_04\3\R01346_R01347_3_var_resampled.tif",
                'R01765':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_09_2\3\R01765_3_var_resampled.tif",
                'R01767':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_09\3\R01767_3_var_resampled.tif"}       
    
      homo_dict = {'R01346':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_04\3\R01346_R01347_3_homo_resampled.tif",
                'R01765':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_09_2\3\R01765_3_homo_resampled.tif",
                'R01767':r"C:\workspace\GLCM\output\new_glcm_rasters\2014_09\3\R01767_3_homo_resampled.tif"}  
                

      rasters = []
      for (key,v), (key1,v1)in zip(var_dict.items(), homo_dict.items())[0:1]:         
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
          
          outFile = r'C:\workspace\GLCM\new_output\gmm_rasters' + os.sep + str(key) + '_GMM_3class_raster.tif'
          rasters.append(outFile)
          CreateRaster(sed_class,gt,outFile)
     
      for (key,v), (key1,v1)in zip(var_dict.items(), homo_dict.items())[0:1]:         
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
               sed_class[ind[k]] = g.predict_proba(vec1[ind[k]].reshape(1, -1))#[:,g.predict(vec1[ind[k]].reshape(1, -1))[0]]
           
          sand = sed_class[:,0] 
           
          gravel = np.add(sed_class[:,1],sed_class[:,2])
          boulders = sed_class[:,3] 
          sand, gravel, boulders =[np.reshape(i,np.shape(read_raster(homo_raster)[0]))for i in [sand, gravel, boulders]]
          #sed_class = np.reshape(sed_class,np.shape(read_raster(homo_raster)[0]))
          
          outFile = r'C:\workspace\GLCM\new_output\gmm_rasters' + os.sep + str(key) + '_sand_GMM3class_proba_raster.tif'
          CreateRaster(sand,gt,outFile)
          outFile = r'C:\workspace\GLCM\new_output\gmm_rasters' + os.sep + str(key) + '_gravel_GMM3class_proba_raster.tif'
          CreateRaster(gravel,gt,outFile)          
          outFile = r'C:\workspace\GLCM\new_output\gmm_rasters' + os.sep + str(key) + '_boulders_GMM3class_proba_raster.tif'
          CreateRaster(boulders,gt,outFile)