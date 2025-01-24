# Changelog

Notable changes to this the AGAGE dataset will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)

## [Unreleased]

## [20250123] - 2025-01-23

This version represents several major changes in the way that AGAGE data are archived and formatted. Data will now be released in a Climate and Forecasting (CF) convention-compliant netCDF format. The archive contains several sets of files:
- a "recommended" file for each site and species, where the relevant ALE/GAGE/AGAGE measurements from a range of instruments have been combined together, with instrument change-over dates specified by the data owner
- "individual instrument" files for each site, species and instruments. This represents the entire public ALE/GAGE/AGAGE dataset, including periods where species have been simultaneously measured in different instruments
- monthly baseline mean mole fractions, representative of monthly averages at a particular site, where conditions have been identified as baseline
- baseline flags as estimated by the Georgia Tech. statistical baseline algorithm

Efforts have been made to restore ALE and GAGE observations and combine them with the AGAGE record in a consistent format, and some early GCMS measurements have also been restored to the archive. 

### Added

- Repeatabilities have been included for ALE and GAGE measurements, based on values inferred from early AGAGE literature (see code repository notes for further details)
- Early GCMS-ADS "Magnum" data have been restored to the public archive
- A "recommended" record for each site and species is provided, so that users do not have to decide which instrument(s) to use for a given time period

### Fixed

- ALE and GAGE timestamps have been converted to UTC (previously local time in some versions)

### Changed

- Data generally released through 2023 (with some exceptions as noted in release schedule in code repository)
- Files are now presented only in netCDF format
- File naming convention has been established:
  - For recommended files: ```agage_<site>_<species>_<version>.nc```
  - For individual instrument files: ```agage-<instrument>_<site>_<species>_<version>.nc```
- Baseline flags are provided in separate files to the observations

### Removed

- Files are not currently provided in ASCII format