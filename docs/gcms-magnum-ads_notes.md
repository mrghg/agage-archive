# Methodology notes
*Matt Rigby, University of Bristol*

with contributions from:

*Ray Wang, GA Tech*

*Joe Pitt, University of Bristol*

---

The first AGAGE automated gas chromatograph-mass spectrometer (GC-MS), based on a Finnigan Magnum Iron Trap coupled to a custom-built (Bristol University) adsorption-desorption system (ADS), was installed at the Mace Head, Ireland (MHD) AGAGE station on October 1994. Here, this system is referred to simply as the "GCMS-Magnum".

In 1997, two new identical ADS GC-MS instruments incorporating Hewlett Packard 5973 quadrupole mass spectrometers were installed at the Mace Head (October 1997) and Cape Grim (November 1997) stations. These instruments superceded the GCMS-Magnum record.

The GCMS-Magnum dataset, and the early MHD HP GCMS-ADS measurements, have not been incorporated into the GCWerks-based data processing workflow that is the basis for the remaining GCMS-ADS and GCMS-Medusa data in this archive. Therefore, we are prepending these data based on a frozen version archived by Georgia Tech. These files are stored in agage/data/data-gcms-magnum.tar.gz and are on the SIO-05 scales.

Repeatabilities were not stored in the archived Magnum files. To estimate repeatabilities, we used the following:
- For compounds where later GCMS-ADS measurments are avaiable at MHD, the mean percentage repeatability was calculated and applied to the Magnum and early ADS data (H-1211, HCFC-141b, HCFC-142b, HFC-152a, HFC-134a, CHCl3, CH2Cl2, CH3Br)
- For CFC-11 and CFC-12, mean percentage repeatabilities at other AGAGE sites that had an ADS were calculated, and the largest mean value was used

These repeatabilities are stored in the file data/agage/gcms-magnum_species.json

