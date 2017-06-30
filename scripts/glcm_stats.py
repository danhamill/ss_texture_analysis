# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 15:32:37 2017

@author: dan
"""
from rasterstats import zonal_stats
import ogr
import pandas as pd
import os
import numpy as np
import platform


def get_subs(shp):
    '''
    Function to get shape file attributes
    Returns: List of substrate attributes
    '''
    ds = ogr.Open(shp)
    lyr = ds.GetLayer(0)
    a=[]
    for row in lyr:
        a.append(row.substrate)
    lyr.ResetReading()
    del ds
    return a

def make_df2(x,metric):
    '''
    Funciton to create empty data frme with specified column name
    Returns: Empty Dataframe
    '''
    df = pd.DataFrame(x,columns=[metric])
    return df

def agg_distributions(stats,in_shp,metric):
    '''
    Funciton to aggregrate the distributions underlying shapefile features
    Returns: Dataframes of the aggregraded distribtuions
    '''
    a = get_subs(in_shp)

    s, g, b = [],[],[]
    n = 0
    for item in stats:
        raster_array = item['mini_raster_array'].compressed()
        substrate = a[n]
        if substrate=='sand':
            s.extend(list(raster_array))
        if substrate=='gravel':
            g.extend(list(raster_array))
        if substrate=='boulders':
            b.extend(list(raster_array))
        n+=1
    del raster_array, substrate, n, item, 

    s_df = make_df2(s,metric)
    g_df = make_df2(g,metric)
    r_df = make_df2(b,metric)
    del s,  g,  b
    return s_df,  g_df, r_df,a


def glcm_summary_stats(k,v,raster,variable,a):
    oName = out_root + os.sep + k +'_' + variable + "_zonalstats" + ".csv" 
    glcm_stats = zonal_stats(v, raster, stats=['min','mean','max','median','std','count','percentile_25','percentile_50','percentile_75'])
    glcm_df = pd.DataFrame(glcm_stats)
    glcm_df['substrate'] = a
    glcm_df.to_csv(oName,sep=',',index=False)
        
if __name__ == '__main__':
    
    
    #Set up root directories for input files and output files
    if platform.system() == 'Windows':
        clone_root = r'c:\workspace\ss_texture_analysis'
        out_root = clone_root + os.sep + "glcm_stats"
    
    #input shapefiles
    shp_dict = {'R01346':clone_root + os.sep + 'shapefiles' + os.sep +'R01346.shp',
                'R01765':clone_root + os.sep + 'shapefiles' + os.sep +'R01765.shp',
                'R01767':clone_root + os.sep + 'shapefiles' + os.sep +'R01767.shp'}  
    fnames = []
    for (k,v) in shp_dict.items()[0:1]:
        
        #Path to textur feature rasters
        ent_file = clone_root + os.sep + 'glcm_rasters' + os.sep + k + '_entropy_resampled.tif'
        var_file = clone_root + os.sep + 'glcm_rasters' + os.sep + k + '_var_resampled.tif'
        hom_file = clone_root + os.sep + 'glcm_rasters' + os.sep + k + '_homo_resampled.tif'
        
        ########################################################################################
        ####            Aggragreate Distibutions for GMM
        ########################################################################################
        #Get mini rasters
        ent_stats = zonal_stats(v, ent_file, stats=['count','mean'], raster_out=True)
        var_stats = zonal_stats(v, var_file, stats=['count','mean'], raster_out=True)
        homo_stats = zonal_stats(v, hom_file, stats=['count','mean'], raster_out=True)
        
        #Aggregrate based on sediment type
        s_ent,  g_ent, r_ent, a = agg_distributions(ent_stats, v,'Entropy')
        s_var,  g_var, r_var = agg_distributions(var_stats, v,'Variance')[0:3]
        s_homo,  g_homo, r_homo = agg_distributions(homo_stats, v,'Homogeneity')[0:3]
        del ent_stats, var_stats, homo_stats
        
        #Merge GLCM properties into sediment type dataframes
        s_df = pd.concat([s_ent,pd.concat([s_var,s_homo],axis=1)],axis=1)
        g_df = pd.concat([g_ent,pd.concat([g_var,g_homo],axis=1)],axis=1)
        r_df = pd.concat([r_ent,pd.concat([r_var,r_homo],axis=1)],axis=1)
        del s_ent, g_ent, r_ent, s_var, g_var, r_var, s_homo, g_homo, r_homo
        
        #Assign numberic sediment type code
        s_df['sedclass'] = 1
        g_df['sedclass'] = 2
        r_df['sedclass'] = 3
        
        #Write distibutions to file
        agg_dist = pd.concat([s_df,pd.concat([g_df,r_df])])
        oName = out_root + os.sep + k + "_aggregraded_distributions.csv"
        fnames.append(oName)
        agg_dist.to_csv(oName,sep=',',index=False)
        
        ########################################################################################
        ####            Calculate Summary Statistics for LSQ
        ########################################################################################
        
        glcm_summary_stats(k,v,ent_file,'ent',a)
        glcm_summary_stats(k,v,hom_file,'hom',a)
        glcm_summary_stats(k,v,var_file,'var',a)
        
        #######NEed to come up with way to get file names for future merge

        
    del k, v, oName, agg_dist, s_df, g_df, r_df,a

    df1 = pd.read_csv(fnames[0],sep=',')  
    df2 = pd.read_csv(fnames[1],sep=',')  
    df3 = pd.read_csv(fnames[2],sep=',')
    
    merge_dist = pd.concat([df1,pd.concat([df2,df3])])

    oName = out_root + os.sep + "merged_aggregraded_distributions_no_outlier.csv.csv"
    merge_dist = merge_dist[(merge_dist['Entropy']>2.5) & (merge_dist['Variance']<12)&(merge_dist['Entropy']<4.6) ]
    merge_dist.to_csv(oName,sep=',',index=False)
    del df1, df2, df3, oName, fnames, merge_dist

    