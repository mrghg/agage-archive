# agage-archive
Code for producing AGAGE archival files

## Installation
Clone this repository using the address in the green ```<> code``` dropdown at the top of the Github page using:

```git clone ...```

Make sure that you have all of the required packages installed. Use the ```requirements.txt``` as a guide. It is recommended that you use conda and conda-forge:

```conda install -c conda-forge <package_name(s)>```

Make sure that this repository is in your path, for example by installing ```conda-build``` and running *from the repo root directory*:

```conda develop .```

## Configuration
Before beginning, run the configuration setup script, which can be accessed as:

```
python agage_archive/util.py
```

Input a descriptive user name when prompted. This function will create a file ```agage_archive/config.yaml```, which contains the default input and output data paths. Note that these paths will be relative to the data/network directory within this repository.

## Usage

We use git and [dvc](https://dvc.org) to track code and data files for this project. Read the [workflow](workflow.md) document for usage details.

## Methodology and references

See the [notes](notes.md) document for methodological details and citations.
