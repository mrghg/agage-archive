# agage-archive
Code for producing AGAGE archival files

## Installation
Clone this repository using the address in the green ```<> code``` dropdown at the top of the Github page using:

```git clone <address>```

It will be easiest to create a new conda environment, ```conda conda create --name agage_env``` then intall ```pip``` using ```conda install pip```.
Alternatively, if avoiding conda, install a virtual environment using ```python -m venv agage_env``` and activate this using ```source agage_env/bin/activate```.

Make sure that you have installed the required dependencies (see ```requirements.txt```), which can be done using ```pip install -r requirements.txt```.

Allow the package to be callable using ```pip install --no-build-isolation --no-deps -e ```. 
If using conda, you can use ```conda develop ```, but note that this functionality is now depricated.

## Configuration
Before beginning, run the configuration setup script, which can be accessed as:

```python agage_archive/config.py```

Input a descriptive user name when prompted. This function will create a file ```agage_archive/config.yaml```, which contains the default input and output data paths. Note that these paths will be relative to the data/network directory within this repository.

## Usage

We use git and [dvc](https://dvc.org) to track code and data files for this project. Read the [workflow](docs/workflow.md) document for usage details.

## Methodology and references

See the [ALE and GAGE notes](docs/ale_gage_notes.md) document for details on how we've processed data from ALE and GAGE, the previous incarnations of AGAGE.

## Reusing this code for processing other datasets

The core functionality behind this code has been designed to be used in other projects. Take a look at the [template-archive](https://github.com/mrghg/template-archive) repository. Note also that processing of flask data requires some specific steps as outlined in the [flasks](docs/flasks.md) document.
