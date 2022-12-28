# <img src="./man/figures/Tropini_Lab_logo.png" align="right" height="150" /> PUPpy - A pipeline to design taxon-specific primers for any defined bacterial community.

### By: Hans Ghezzi - Tropini Lab

# Project 

PUPpy (Phylogenetically Unique Primers in python) is a computational pipeline developed to design taxon-specific primers within a defined bacterial community. PUPpy can design both strain-specific primers, which selectively amplify each member of the community, and group-specific primers, which selectively amplify user-selected members. Primers designed with the pipeline can be used to assess the presence/absence of bacteria in samples through PCR, as well as quantify their abundance via qPCR. 

# Installation

PUPpy can be installed thorugh Conda for mac and Linux, or by installing the required dependencies.

```sh 
# Install from Conda
<ADD CODE HERE>
```

If this were to fail, you can create the required Conda environment from this GitHub directory.

```sh
# Deactivate Conda environment
conda deactivate

# Clone PUPpy GitHub directory
git clone https://github.com/Tropini-lab/PUPpy_pipeline.git

# Create conda environment
cd PUPpy_pipeline
conda env create -f environment.yml
conda activate PUPpy_pipeline
```

Alternatively, you can manually set up a conda environment with the individual dependencies:

```sh
# 1) Create a Conda environment
conda create -n PUPpy_pipeline python=3.10.6

# 2) Activate the Conda environment
conda activate PUPpy_pipeline

# 3) Configure Conda channels
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# 4) Install dependencies:
# Conda
## MMseqs2 v14.7e284
conda install -c "bioconda/label/cf201901" mmseqs2
## Pandas 1.5
conda install pandas=1.5
## BioPython v1.80
conda install -c "conda-forge/label/cf201901" biopython
## Dask v0.15.2
conda install -c "conda-forge/label/cf201901" dask
## R tidyverse v1.3.1
conda install -c r r-tidyverse
## R readr v2.1.3
conda install -c "conda-forge/label/cf201901" r-readr
## R stringi v1.7.8
conda install -c "conda-forge/label/cf201901" r-stringi
# Pip
## primer3-py 0.6.1
pip install primer3-py
## colorama 0.4.6
pip install colorama
## psutil 5.9.4
pip install psutil
```

# Simple usage

PUPpy consists of 2 main steps: 1) aligning the input genes and 2) designing taxon-specific primers.

The alignment step must always be run for any new defined bacterial community:

```python
python alignments.py -c input_CDS_dir -o OUT_dir
```

The second step is where users choose whether to design taxon-specific primers unique to individual members or shared by groups of the bacterial community.

```python
python primerDesign.py -t Target_Species_CDSes -i input_alignments.tsv -o OUT_primerDir
```
