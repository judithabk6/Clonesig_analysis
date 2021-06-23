# Clonesig_analysis
This repository contains all the scripts used to provide the results for the CloneSig article.

Abécassis, Judith, Fabien Reyal, and Jean-Philippe Vert. "CloneSig: Joint inference of intra-tumor heterogeneity and signature deconvolution in tumor bulk sequencing data." BioRxiv (2019): [825778](https://www.biorxiv.org/content/10.1101/825778v2).


The folder ```signature_code``` contains most of the code, with the master script in the bash script ```run_all.sh```. It is not advised to execute it all at once, considering the run time, and the paths to change (python virtual environments, paths to executables etc).

The folder ```notebooks``` contains notebooks and scripts in which the final analyses are generated. The tables used to generate those results are in the folder ```result_tables```.

The folder ```external_data``` contains data from external sources that were useful to carry out the project. A specific readme is included to provide the origin of each file.

It should be noted that more complete and up-to-date instructions of the CloneSig python package can be found [here](https://github.com/judithabk6/clonesig). Feel free to get in touch or to open an issue if you have any question or suggestion.


#### Installation and dependencies
this project uses many packages. Here is a list of the most important ones. Requirement (```python27_requirements.txt``` and ```python36_requirements.txt```) files are provided to use with pip and setup two virtual environnments for Python 2.7 and 3.6 (or conda for some specific packages, in particular ```mkl``` and ```PyClone```). Here is a list of the main ones, including some that can be installed with pip (specified). Code was run with R version 3.3.2, and Python 3.6.8 (Anaconda installation), except to run PyClone (Python 2.7.9). All computation was performed under a Centos distribution, with torque as scheduler.

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
- [TrackSigFreq 23b2f3f](https://github.com/morrislab/TrackSigFreq/tree/23b2f3f75b344d18d1df6817f3492f8b80047500)
- [DPClust 75f5d7e](https://github.com/Wedge-lab/dpclust/tree/75f5d7ef1e3e53585f86801fde76dd4c4aa86324)
- [PhylogicNDT c229cec](https://github.com/broadinstitute/PhylogicNDT/tree/c229cec570169b6e710e9157c9e102ce37a454cd)


##### statistical packages
- [lifelines](https://lifelines.readthedocs.io/en/latest/Quickstart.html)

All scripts should be run from the root of the folder. Total run time is several weeks on a cluster with 60 CPUs.


