# Optical data processing notes
*Matt Rigby, University of Bristol*

with contributions from:

*Joe Pitt, University of Bristol*

---

Optical instruments have been processed and combined with the earlier ALE/GAGE/AGAGE records for some species (e.g., methane). Here, optical instruments refer to laser-based systems such as Picarro and LGR systems that are based on Cavity Ring-down Spectrometry (CRDS). 

## Naming

Several AGAGE sites have now used multiple Picarro instruments. In the GCWerks files, these are output with filenames "...Picarro-1...", "...Picarro-2...", etc. In this code, we keep those names and treat them as individual instruments (c.f., GCMD, GCMS-Medusa). Therefore, if multiple instruments are to be combined, they should be in separate files/columns in the release schedules and data combination files with the instrument names "Picarro-1", "Picarro-2", etc.

## Resampling

For sites with only one inlet, or where inlets have changed only occasionally, the high-frequency optical data have been resampled to hourly. However, for sites like Tacolneston (TAC), the inlet is switched between three heights every 20 minutes, and therefore hourly resampling would lose information from some of the heights, or average over heights. To address this, data are grouped into ~20 minute bins, reflecting the sampling period for each inlet. Whether to resample or group be inlet height is decided by the inlet switching frequency in the function ```convert.resample```. Resampling is done by the ```resampler``` function, and grouping by the ```grouper``` function. The latter is computationally expensive, so processing is slow when this function is used.