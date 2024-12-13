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
snakemake environment
```

To create the tcr-antigen-prediction mamba environment.

If **only mamba is installed**, manually run the following steps to install the project dependencies (including Snakemake):

```
mamba env create -f environment.yml
mamba run -n tcr-antigen-prediction python -m pip install .

git submodules init
git submodule update --recursive
mamba_prefix=$(mamba run -n tcr-antigen-prediction mamba info --json | jq '."env location"' | sed s/\"//g)
python_version=$(mamba run -n tcr-antigen-prediction python --version | cut -d " " -f2 | cut -d "." -f1-2)
cp -r third_party/anarci $mamba_prefix/lib/python$python_version/site-packages
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
snakemake test
```
