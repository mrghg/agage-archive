# Updating the AGAGE public archive

This repository uses two types of version control:
- git is used to track versions of the code
- [dvc](https://dvc.org) is used to track versions of the processed data

Make sure that you have git installed on your system before you start.

The dvc and dvc-gcloud (Google Cloud plugin) are listed in the requirements file, and can be installed from conda, using the conda-forge channel.

## Initial setup

It will probably be easiest to create a new conda environment, ```conda conda create --name agage_env``` then install ```pip``` using ```conda install pip```.
Alternatively, if avoiding conda, install a virtual environment using ```python -m venv agage_env``` and activate this using ```source agage_env/bin/activate```.

Make sure that you have installed the required dependencies (see ```requirements.txt```), which can be done using ```pip install -r requirements.txt```.

Allow the package to be callable using ```pip install --no-build-isolation --no-deps -e . ``` and run ```python -c 'from agage_archive import util; util.setup()'```, to set the desired input and output paths.

### DVC initial setup

Firstly, check that DVC is working by running from the terminal:

```dvc status```

It should say something like:

```
data/agage/agage-public-archive.zip.dvc:                              
        changed outs:
                not in cache:       data/agage/agage-public-archive.zip
```

The first time you use DVC, you will need to provide some credentials for the remote file hosting. We are using Google Drive for this. We are using a Google Drive service account.  **These credentials must not be shared outside of the AGAGE team. Please ensure that the authentication file is NOT committed to the repository or shared online.**

You will need to contact Matt Rigby for the authentication json file. Put this file somewhere on your local machine. You can put it in the repository folder if you like, where it will be ignored by the version control (assuming that you don't change the name).

Run the following command from a terminal in the repository directory, replacing ```/path/to``` with the path to the location of the json file:

```dvc remote modify myremote --local gdrive_service_account_json_file_path /path/to/agage-gdrive.json```

Note that the ```--local``` flag is extremely important here. It adds a file ```.dvc/config.local``` containing the location of the authentication file. This local config file is excluded from the repository.

**If the ```agage-gdrive.json``` file is accidentally pushed to the repository or shared outside of the AGAGE team, please contact Matt Rigby immediately.**

You will find out whether the authentication has worked the first time you do ```dvc pull```.

### Note to admin

Authentication settings are through the [Google Drive API](https://console.developers.google.com/). See instructions on authentication and service accounts on the [DVC gdrive docs](https://dvc.org/doc/user-guide/data-management/remote-storage/google-drive).

Need to set:

```dvc remote modify myremote gdrive_use_service_account true```

Which should then be in .dvc/config

## Usage

Firstly, make sure that the git repository is up to date (make sure you know how to use git, but assuming you do, just do ```git pull``` to get the latest version).

You should now be able to run:

```dvc pull```

to get the latest version of the public archive (file ```data/agage/agage-public-archive.zip```).

If you want to update the archive:

1. Download the latest data from GCWerks, and place it at the locations specified in your ```config.yaml``` file.
2. Modify the data specification files as required, in ```data/agage```:
```
data_release_schedule.xlsx
data_combination.xlsx
data_exclude.xlsx
scale_defaults.csv
```
3. Re-run ```run.py```, or any other processing functions that you'd like to run. 
4. Thoroughly check the processed data files by examining the netcdf files stored in the newly created ```data/agage/agage-public-archive.zip``` (the functions in ```io.py``` can automatically extract files from the zip archive, given the appropriate site, instrument, etc.), and using the ```notebooks/visualise.ipynb``` data viewer.
5. Once you're happy, run: 
```
dvc add data/agage/agage-public-archive.zip
git commit -am "COMMENT DESCRIBING THE CHANGE YOU'VE MADE"
dvc push
git push
```

## Copying relevant input files

You will need to copy over the GCWerks output folders ```data-nc.zip``` and ```data-gcms-nc.zip``` from whereever you pull the data to your local ```data/archive``` folder (and zip the two folders if not zipped already).