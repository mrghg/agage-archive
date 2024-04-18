
#  AGAGE Data README file
This archive contains the complete AGAGE/ALE/GAGE dataset. 

The archive is compiled from early ALE/GAGE files combined with AGAGE files output from the GCWerks software. The code for processing and standardising these data is at https://github.com/mrghg/agage-archive. The ALE/GAGE files were created by Derek Cunnold's team at Georgia Tech (GA Tech.) in the early 1990s, and are archived at https://github.com/mrghg/agage-archive/tree/main/data/agage/ale_gage_sio1993. These data were themselves converted from earlier scales to the Scripps Institution of Oceanography (SIO) 1993 scale (SIO-93; see Cunnold et al., 1994 and Prinn et al., 2000). The GA Tech. files are reprocessed here and converted to more recent scales. See the Methodology notes in the agage-archive repository for details.

The AGAGE data come from netCDF files generated by GCWerks, with some minor modifications using the agage-archive code.

The code for reprocessing the AGAGE data is currently maintained by the University of Bristol. Contact matt.rigby@bristol.ac.uk with errors and questions. For questions regarding individual data sources, please contact the relevant AGAGE station PIs (contact details are contained in file metadata).

# Format and Structure

This archive consists of a set of netCDF files. For information on the netCDF file format see: https://www.unidata.ucar.edu/software/netcdf/

Filenames follow the convention:

```network-instrument_sitecode_species_version.nc```

Below is an example of the archive structure showing only methane data at two sites, CGO and ZEP:

```
├── README.md
├── event
│   └── ch4
│       ├── AGAGE-PICARRO_ZEP_ch4_20240208v1.nc
│       ├── AGAGE-combined_CGO_ch4_20240208v1.nc
│       ├── baseline_flags
│       │   ├── AGAGE-PICARRO_ZEP_ch4-git-baseline_20240208v1.nc
│       │   └── AGAGE-combined_CGO_ch4-git-baseline_20240208v1.nc
│       └── individual
│           ├── AGAGE-GAGE-GCMD_CGO_ch4_20240208v1.nc
│           ├── AGAGE-GCMD_CGO_ch4_20240208v1.nc
│           ├── AGAGE-PICARRO_CGO_ch4_20240208v1.nc
│           └── baseline_flags
│               ├── AGAGE-GAGE-GCMD_CGO_ch4-git-baseline_20240208v1.nc
│               ├── AGAGE-GCMD_CGO_ch4-git-baseline_20240208v1.nc
│               └── AGAGE-PICARRO_CGO_ch4-git-baseline_20240208v1.nc
└── monthly
    └── ch4
        ├── AGAGE-PICARRO_ZEP_ch4-monthly_20240208v1.nc
        ├── AGAGE-combined_CGO_ch4-monthly_20240208v1.nc
        └── individual
            ├── AGAGE-GAGE-GCMD_CGO_ch4-monthly_20240208v1.nc
            ├── AGAGE-GCMD_CGO_ch4-monthly_20240208v1.nc
            └── AGAGE-PICARRO_CGO_ch4-monthly_20240208v1.nc
```

At the first level, the archive is organised into high-frequency (~hourly) files, in a folder called "event", and monthly baseline averages ("monthly"). For the latter, baselines have been estimated using the AGAGE statistical pollution algorithm (see O'Doherty et al., 2001).

Within each of these top-level folders are directories containing the data for each species.

Within the species directory, if multiple instruments have measured a species at a single site, there is a **'combined'** .nc file, where individual instruments have been combined to form a continuous record (e.g., for CGO in the above example). Choices on when to switch instruments are provided in the data_combination.xlsx spreadsheet in the agage-archive repository.

Where data have been combined, the individual data files for each instrument can be found in the sub-directory ```individual```.

If a species is only measured on one instrument (e.g., for ZEP in the above), the indiviudal instrument files are in the species folder and do not appear in the ```individual``` folder. Therefore, most users will only be interested in the files contained at the top level of the species file. The "individual" folder is required only if users are interested in indiviudal instrument performance and intercomparison.

The ```baseline``` sub-folders contains the baseline flags associated with each type of ```event``` file.

# AGAGE Data Statement:

These data are made available to the scientific community in the belief that their wide use will lead to new scientific insights. The availability of these data does not constitute publication of the data. AGAGE relies on the ethics and integrity of the user to ensure that the AGAGE scientists receive fair credit for their work. If the data are obtained for potential use in a publication or presentation, AGAGE should be informed at the outset of this work. If the AGAGE  data are essential to the work, or if an important result or conclusion depends on the AGAGE data, co-authorship may be appropriate. This should be discussed at an early stage with the AGAGE contacts listed below. Manuscripts using the AGAGE data should be sent to the AGAGE contacts for review before they are submitted for publication so we can ensure that the quality and limitations of the data are accurately represented. Every effort is made to produce the most accurate and precise measurements possible. However, we reserve the right to make corrections to the data based on recalibration of standard gases or for other reasons deemed scientifically justified. We are not responsible for results and conclusions based on use of these data without regard to this warning.

# AGAGE Data Reciprocity Agreement

Use of these data implies an agreement to reciprocate. Laboratories making similar measurements agree to make their own data available to the general public and to the scientific community in an equally complete and easily accessible form. Scientists  are encouraged to make available to the community, upon request, their own modelling tools used in the interpretation of the AGAGE data, namely well documented model code, transport fields, and additional information necessary for other scientists to repeat the work and to run modified versions.

# Acknowledgements

Please thank the listed contacts (who are not already co-authors) for each of the stations whose data is used in your paper. Please state: “AGAGE is supported by NASA (USA) grants to MIT and SIO, UK government and NOAA (USA) grants to Bristol University; CSIRO and BoM (Australia); FOEN grants to Empa (Switzerland); NILU (Norway); SNU (Korea); CMA (China); NIES (Japan); and Urbino University (Italy)”.

# Citation of AGAGE Data
 
(1) Obtain the DOI of the specific version of the AGAGE data that you have used.

(2) R.G. Prinn, R.F. Weiss, J. Arduini, T. Arnold, H.L. DeWitt, P.J. Fraser, A.L. Ganesan, J. Gasore, C.M. Harth, O. Hermansen, J. Kim, P.B. Krummel, S. Li, 
Z. M. Loh, C.R. Lunder, M. Maione, A.J. Manning, B.R. Miller, B. Mitrevski, J. Mühle, S. O’Doherty, S. Park, S. Reimann, M. Rigby, T. Saito, P.K. Salameh, 
R. Schmidt, P.G.  Simmonds, L.P. Steele, M.K. Vollmer, R.H. Wang, B. Yao, Y. Yokouchi, D. Young, and L. Zhou: History of chemically and radiatively important 
atmospheric gases from the Advanced Global Atmospheric Gases Experiment (AGAGE), Earth Syst. Sci. Data, 10, 985-1018, 2018, https://doi.org/10.5194/essd-10-985-2018. 

(3) Most recent AGAGE paper publishing the relevant data (e.g., for a specific species), as listed on AGAGE Publications.

# Licence

Distributed under CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/deed.en)