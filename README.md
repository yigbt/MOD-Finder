# MOD-Finder

**MOD-Finder** (Multi-Omics Data set Finder) is an easy-to-use R Shiny
app to efficiently search for compound-related data sets on the
transcriptomic, proteomic, and metabolomic layer.

The webapp is available at [UFZ](https://webapp.ufz.de/mod_finder/),
but the code can be run as a local stand-alone application. To do so,
just clone this repository and follow the instructions described under
'Prerequisites'. 




## Prerequisites

Before the application can be started locally, some requirements have
to be met. First, **MOD-Finder** was developed in a specific conda
environment to avoid changing R package to infer with the
functionality and stability of **MOD-Finder**. Second, **MOD-Finder**
relies on public online sources that are subject to changes once in a
while. But since not all of these sources provide sufficient access
via RESTful APIs for example, **MOD-Finder** needs several local
database files that have to be generated in advance.


### Conda environment

To ensure reproducibility and functionality of **MOD-Finder**, we
encourage all users to run the app in the provided conda environment.
The environment can be created and activate by

```bash
conda env create -f mod_finder_conda.yml
source activate mod_finder
```

### Creating local database files

These database files are essential for **MOD-Finder** in order to function properly.
Before the apps is started for the first time, three different database files **have to be present**:

* comptox.RData
* CTD_chemicals.RData
* GEOmetadb.sqlite

The first file `comptox.RData` is essential for the *digital*
identification of chemicals of interest and contains an ID mapping of
more than ??? thousand chemicals, including for example _Pubchem IDs_,
_CAS RN_, and _Comptox IDs_. `CTD_chemicals` is the preprocessed
source for retrieving chemical-specific information about effects on
gene expression, pathway perturbations, and associations to diseases.
`GEOmetadb.sqlite` contains the current version of the NCBI Geo
database in an SQL format. 

To create these database files in advance type the following commands:

```bash
cd /path/to/app/
mkdir data/
bash scripts/create_local_databases.sh data/
```


## How to Cite

When you find **MOD-Finder** useful, please make sure to cite us accordingly:
_still to come_
