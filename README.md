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

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) are used to setup and manage the dependencies and workflows in this project. At a minimum, Conda is required to install the dependencies for the project (a version of Snakemake will be installed as a dependency).

If **both Conda and Snakemake are installed**, simply run:

```
snakemake environment
```

To create the tcr-antigen-prediction conda environment.

If **only conda is installed**, manually run the following steps to install the project dependencies (including Snakemake):

```
conda env create -f environment.yml
conda run -n tcr-antigen-prediction python -m pip install .

git submodules init
git submodule update --recursive
conda_prefix=$(conda run -n tcr-antigen-prediction conda info --json | jq .default_prefix | sed s/\"//g)
python_version=$(conda run -n tcr-antigen-prediction python --version | cut -d " " -f2 | cut -d "." -f1-2)
cp -r third_party/anarci $conda_prefix/lib/python$python_version/site-packages
```

The environment can then be activated:

```
conda activate tcr-antigen-prediction
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