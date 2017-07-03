####

Welcome to the repository for GLCM based texture analsis of side scan sonar echograms collected with a recreational-grade system.  

#### Contents
This repostitory contains all of the data and scripts used to prepare the figures in _Alluvial mappint by automated texture segmentation of recreational-grade side scan sonar imagery_ submitted to Environmental Modelling and Sofware. This repository contains a collection of files and scripts required to peform the analysis.  

#### Organization

This repository is organized as follows:
* `/ss_rasters/` contains georeferenced side scan sonar echograms
* `/shapefiles/` contains shapefiles of the vidisually identified sediment patches
* `/sedclass_rasters/` contains georeferenced sediment classification maps
* `/scripts/` contains all of the scripts requried to create sediment classication maps
* `/glcm_stats/` contains csv files of GLCM distribuions and GLCM summary statistics
* `/glcm_rasters/` contains georeferenced GLCM texture feature rasters

#### Dependencies
All of the scripts were developed using python 2.7.11 in a windows 10 enviroment.  I used the Anaconda distribution 4.0.0 (64 bit) with the MSC v.1500 64 bit (AMD64) compiler.  The following dependencies are required:


* [gdal 1.11.4 python bindings](http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal)
* [gdal 1.11.4](http://www.gisinternals.com/release.php)
* [rasterstats](https://github.com/perrygeo/python-rasterstats)
* numpy
* pandas
* sklearn
* [scikit-image](https://github.com/scikit-image/scikit-image)
* [joblib](https://github.com/joblib/joblib)
* [pyproj 1.9.4](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyproj)


#### Workflow
All of the continious side scan sonar recordings were processed using [PyHum](https://github.com/dbuscombe-usgs/PyHum).  In the interest of space, I have not included any of the binary side scan sonar files, intermediate PyHum files, or georeferenced point clouds.  If any of those files are of interest, please contact me and I will provided them outside of this repository. 

To start, I reccomentd cloning this repository to `c:\workspace`.  If you want to clone the repository to a differnet directory, there is a variable `clone_root` at the begining of each python script where you can indicate where the appropriate directory.

Begining with the side scan sonar echogram rasters in `/ss_rasters/`, you will first need to calculate GLCM texture features using `/scripts/GLCM_calc.py`.  This script will produce georeferenced GLCM texture features in the directory `/output/glcm_rasters/`. 

```python
python glcm_calc.py
```
Next you will have to use the shapefiles provided in `/shapefiles/` to calculate sediment type statistics.  Aggregraded distibutions and summary statistic CSV file will be save to `/glcm_stats/`.

```python
python glcm_stats.py
```
The statistics saved in `/glcm_stats/` are used to calibrate all of the automated texture segmentaion algorithims.  All of the sediment classification rasters will be output to `/sedclass_rasters/`.  There are individual python scripts for each of the texture segmentation methods in scipts.  

```python
python LSQ.py
python gmm2.py
python gmm4.py
```

