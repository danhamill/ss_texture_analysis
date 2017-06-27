####

Welcome to the repository for GLCM based texture analsis of side scan sonar echograms collected with a recreational-grade system.  

#### Contents
This repostitory contains all of the data and scripts used to prepare the figures in _Coarse-resolution alluvial mappint by automated texture segmentation of recreational-grade side scan sonar imagery_ submitted to Environmental Modelling and Sofware. This repository contains a collection of files and scripts required to peform the analysis.  

#### Organization

This repository is organized as follows 

#### Dependencies
All of the scripts were developed using python 2.7.11 in a windows 10 enviroment.  I used the Anaconda distribution 4.0.0 (64 bit) with the MSC v.1500 64 bit (AMD64) compiler.  The following dependencies are required:


* [gdal 1.11.4](http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal)
* [rasterstats](https://github.com/perrygeo/python-rasterstats)
* numpy
* pandas
* sklearn
* [scikit-image](https://github.com/scikit-image/scikit-image)
* [joblib](https://github.com/joblib/joblib)
* [pyproj 1.9.4](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyproj)


#### Workflow
All of the continious side scan sonar recordings were processed using [PyHum](https://github.com/dbuscombe-usgs/PyHum).  In the interest of space, I have not included any of the binary side scan sonar files, intermediate PyHum files, or georeferenced point clouds.  If any of those files are of interest, please contact me and I will provided them outside of this repository.  However, I have included some auxillary scripts if you need the scripts to resample the point clouds to a reqular grid and convert them to raster format in '/scripts/extra_scripts/'.  Note, you will need install [pyresample 1.1.4](http://pyresample.readthedocs.io/en/latest/). The latest version (1.1.6) at the time of writing this readme does not play nicely with scipy.  Hence, the previous version.

Begining with the side scan sonar echogram rasters in '/ss_rasters/`, you will first need to calculate GLCM texture features using `/scripts/GLCM_calc.py`.  This script will produce georeferenced GLCM texture features in the directory `/output/glcm_rasters/`.  Before any zonal statistics are calculated the GLCM texture features need to be resampled to a 3-meter resolution.  There are many ways to do this, but I reccomend using the command line utlilty [gdalwarp](http://www.gdal.org/gdalwarp.html).

```
gdalwarp -tr 3 3 -srcnodata -99 -dstnodata -99 -r average c:/workspace/ss_texture_analysis/output/glcm_rasters/R01346_R01347_3
