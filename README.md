# TCR Antigen Prediction

Models and methods for predicting TCR antigen specificity.

## Installation and Setup

To install this project, first clone the repo:

```
git clone git@github.com:benjiemc/tcr-antigen-prediction.git
cd tcr-antigen-prediction/

# Or using HTTPS
git clone https://github.com/benjiemc/tcr-antigen-prediction.git
cd tcr-antigen-prediction/
```

> ### Using MODELLER to fix missing residues and atoms
>
> Parts of the data pipeline for this project use [MODELLER](https://salilab.org/modeller/) to fix structures missing atoms and residues. A license is required to use MODELLER. Please [register for a license](https://salilab.org/modeller/registration.html) and then set the license key as an evironment variable using the following command: `export KEY_MODELLER='XXXXXX'`, replacing XXXXXX with your assigned license key.
>
> Alternatively, the use of MODELLER can be manually disabled in the workflow steps.

[Mamba](https://mamba.readthedocs.io/en/latest/index.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) are used to setup and manage the dependencies and workflows in this project. At a minimum, Mamba is required to install the dependencies for the project (a version of Snakemake will be installed as a dependency).

If **both Mamba and Snakemake are installed**, simply run:

```
snakemake environment --cores 1
```

To create the tcr-antigen-prediction mamba environment.

If **only mamba is installed**, manually run the following steps to install the project dependencies (including Snakemake):

```
# Create Environment
git submodule init
mamba env create -f environment.yml
mamba run -n tcr-antigen-prediction python -m pip install .

# Download Stitchr Data
mamba run -n tcr-antigen-prediction stitchrdl -s human
mamba run -n tcr-antigen-prediction stitchrdl -s mouse

# Install ANARCI
git submodule update --remote third_party/anarci
mamba_prefix=$(mamba run -n tcr-antigen-prediction mamba info --json | jq '."env location"' | sed s/\"//g)
python_version=$(mamba run -n tcr-antigen-prediction python --version | cut -d " " -f2 | cut -d "." -f1-2)
cp -r third_party/anarci $mamba_prefix/lib/python$python_version/site-packages

# Install NetTCR-clone
git submodule update --remote third_party/NetTCR
mamba run -n tcr-antigen-prediction pip install third_party/NetTCR
```

The environment can then be activated:

```
mamba activate tcr-antigen-prediction
```

### Installing for Development

To install the required dependencies for development, run the above and additionally run the command:

```
pip install -e '.[develop]'
```

This will install the required development dependencies. To verfiy the installation, run:

```
snakemake test --cores 1
```

## Database Notes

Several sequence databases were manually downloaded to make up the datasets for this project:

1. IEDB
    * URL: https://www.iedb.org/
    * Date: 2025/01/14
    * Selection Criteria:
      * Assay: T cell, MHC Ligand, Outcome - Positive
      * Epitope: Any
      * MHC Restriction: Any
      * Host: Any
      * Disease: Any
      * Reference: Any
    * Export T Cell Receptors with Epitopes (single headers as csv)
    * Result: 224498 sequences downloaded

2. ITRAP
    * URL: https://doi.org/10.11583/DTU.22645342.v1
    * Date: 2025/01/14
    * Complete download of TCR (highest quality) dataset (2833 sequences)

3. McPas-TCR
    * URL: https://friedmanlab.weizmann.ac.il/McPAS-TCR/
    * Date: 2025/01/10
    * Download complete database (40779 sequences)

4. VDJdb
    * URL: https://vdjdb.cdr3.net/
    * Date: 2025/01/10
    * Selection Criteria:
      * CDR3 - Species: Human, Monkey, Mouse
      * CDR3 - Gene (chain): TRA, TRB
      * MHC - Class: MHCI, MHCII
      * Meta - Assay Type: Multimer sorting, Culture-based, Other
      * Meta - Sequencing: Sanger, High-throughput, Single-cell
      * Meta - Spurious CDR3: Include non-canonical, Include unmapped V/J
    * Result: All sequences downloaded (120849 of 120849)
