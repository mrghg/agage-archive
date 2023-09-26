# Updating the AGAGE public archive

This repository uses two types of version control:
- git is used to track versions of the code
- [dvc](https://dvc.org) is used to track versions of the processed data

Make sure that you have git installed on your system before you start.

The dvc and dvc-gcloud (Google Cloud plugin) are listed in the requirements file, and can be installed from conda, using the conda-forge channel.

## Initial setup

Make sure that you have installed the required dependencies (see ```requirements.txt```), and run ```setup.py```, to set the desired input and output paths.

### DVC initial setup

Firstly, check that DVC is working by running from the terminal:

```dvc status```

It should say something like:

```
data/agage/agage-public-archive.zip.dvc:                              
        changed outs:
                not in cache:       data/agage/agage-public-archive.zip
```

The first time you use DVC, you will need to provide some credentials for the remote file hosting. We are using Google Drive for this. You will need to set up a `client_id` and `client_secret` to access the remote storage. You will need to contact Matt Rigby and provide a Google Drive account name to access these. **These credentials must not be shared outside of the AGAGE team.**

Run the following commands from a terminal in the repository directory:

```dvc remote modify --local myremote gdrive_client_id 'CLIENT_ID_FROM_MATT'```

```dvc remote modify --local myremote gdrive_client_secret 'CLIENT_SECRET_FROM_MATT'```

Note that the ```--local``` flag is extremely important here. It adds a file ```.dvc/config.local``` containing the authentication credentials, which is excluded from the git repository. If it is not used, these credentials are added to ```.dvc/config```, which is included in the git repo, and can then potentially be visible to other users. **Please ensure that you do not accidentally commit these details to ```.dvc/config```**.

The first time that you do ```dvc pull``` or similar interaction with the remote repository, you will be taken to a Google drive authentication screen on your browser. Follow the instructions, and you should be granted access.

### Note to admin

User will need to be added as a test user under OAuth Consent Screen on the [Google Drive API](https://console.developers.google.com/) console, **and** they will needed to be granted access to the Google drive folder.

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
