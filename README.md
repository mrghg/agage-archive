# agage-archive
Code for producing AGAGE archival files

## Installation
Clone this repository using the address in the green ```<> code``` dropdown at the top of the Github page using:

```git clone <address>```

It will be easiest to create a new conda environment, ```conda conda create --name agage_env``` then intall ```pip``` using ```conda install pip```.
Alternatively, if avoiding conda, install a virtual environment using ```python -m venv agage_env``` and activate this using ```source agage_env/bin/activate```.

Make sure that you have installed the required dependencies (see ```requirements.txt```), which can be done using ```pip install -r requirements.txt```.

Allow the package to be callable using ```pip install --no-build-isolation --no-deps -e . 

## Configuration
Before beginning, run the configuration setup script, which can be accessed as:

```python agage_archive/util.py```

Input a descriptive user name when prompted. This function will create a file ```agage_archive/config.yaml```, which contains the default input and output data paths. Note that these paths will be relative to the data/network directory within this repository.

## Usage

We use git and [dvc](https://dvc.org) to track code and data files for this project. Read the [workflow](workflow.md) document for usage details.

## Methodology and references

See the [notes](notes.md) document for methodological details and citations.
