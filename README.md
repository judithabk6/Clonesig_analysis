# Clonesig_analysis
This repository contains all the scripts used to provide the results for Clonesig publication.


The folder ```signature_code``` contains most of the code, with the master script in the bash script ```run_all.sh```. It is not advised to execute it all at once, considering the run time, and the paths to change (python virtual environments, paths to executables etc).

The folder ```notebooks``` contains notebooks and scripts in which the final analyses are generated. The tables used to generate those results are in the folder ```result_tables```.

The folder ```external_data``` contains data from external sources that were useful to carry out the project. A specific readme is included to provide the origin of each file.

It should be noted that more complete and up-to-date instructions of the CloneSig python package can be found [here](https://github.com/judithabk6/clonesig). Feel free to get in touch or to open an issue if you have any question or suggestion.


#### Installation and dependencies
this project uses many packages. Here is a list of the most important ones. Requirement files are provided to use with pip (or conda for some specific packages, in particular ```mkl``` and ```PyClone```). Here is a list of the main ones, including some that can be installed with pip (specified). Code was run with R version 3.3.2, and Python 3.6.8 (Anaconda installation), except to run PyClone (Python 2.7.9). All computation was performed under a Centos distribution, with torque as scheduler.

##### data download packages
- [gdc transfer tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)

##### data preprocessing
- [segment_liftover 0.951](https://pypi.org/project/segment-liftover/) (python >= 3.6 - installed by pip)

##### ITH and signature methods
- [PyClone 0.13.0](https://bitbucket.org/aroth85/pyclone/wiki/Installation) (Python 2.7)
- [Sciclone 1.1](https://github.com/genome/sciclone)
- [Tracksig](https://github.com/morrislab/TrackSig)
- [Ccube 205f5e7b89](https://github.com/keyuan/ccube/tree/205f5e7b895fd302019faced3f4acf1a3f15e778)
- [palimpsest](https://github.com/FunGeST/Palimpsest)
- [CloneSig](https://github.com/judithabk6/clonesig)
- [deconstructSigs 9bbaf15387](https://github.com/raerose01/deconstructSigs/tree/9bbaf15387e1a6221b4437523d12dd950eea80e1)


##### statistical packages
- [lifelines](https://lifelines.readthedocs.io/en/latest/Quickstart.html)

All scripts should be run from the root of the folder.
