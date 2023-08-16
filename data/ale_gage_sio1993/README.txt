The files in this folder are the ALE/GAGE data prepared by Derek Cunnold at GA Tech in the early 1990s. They were converted to the SIO-93 scale (CSIRO-94 for CH4) from the original ALE scales. Data are dry air mole fractions. The file format is described in notes.md, in the root director of this repository.

There are numerous issues with these files, which are discussed in the "notes.md" file included in the root directory of this repository. Therefore, these files should not be used in isolation from the various processing routines and other data sources compiled in this repository.

The files can be read and cleaned using the function io.read_ale_gage.