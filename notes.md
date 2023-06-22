# Methodology notes
*Matt Rigby, University of Bristol*

with contributions from:

*Ray Wang, GA Tech*

*Paul Krummel, CSIRO*

---

Some features of the ALE/GAGE data have to be reconstructed from literature or personal communication. Here are some bits of information that have been provided to me, or dug out of key papers at the time.

## File format
The ALE/GAGE files were created by Derek Cunnold and his team at Georgia Tech (GA Tech.) in the early 1990s. These data were all converted from earlier scales to the Scripps Institution of Oceanography (SIO) 1993 scale (SIO-93; see Cunnold et al., 1994 and Prinn et al., 2000).

Data are in *local time* at each site. One file per site, per month.

Columns refer to:

- 1st column, day of month
- 2nd column, time in hhmm (hour, minute)
- 3rd column, "ABSDA", days elapsed since the start of ALE program
- 4th through last column is the dry air mole fraction of measured compounds. -99.0 indicates missing values. The "P" after mole fraction indicates pollution event according to the GA Tech statistcal pollution algorithm.

F-11S and F-11P refer to silicone or Porasil chromatographic columns. The "S" column is preferred.

## Scale conversions

|Species|SIO-98/SIO-93|SIO-05/SIO-98|TU-87/CSIRO-94|
|--|--|--|--|
|CFC-11|1.0082|0.9945||
|CFC-12|1.0053|1.0000||
|CFC-113|1.0038|1.0042||
|CH3CCl3|1.0184|0.9957||
|CCl4|1.0090|1.0000||
|CH4*|||1.0119|
|N2O**|1.0058|1.0000||


*CH4 note: the factor of 1.01069 was originally used to adjust CH4 to Japanese (TU-1987) scale (when we converted other compounds to SIO-98 scale), it was changed to 1.0119 after Paul Steel (CSIRO) updated the factor (when we updated other species to SIO-2005 scale).
 
**N2O, in addition to applying the factor from SIO-93 to SIO-98, the value also needs to be divided by a time-dependent polynomial function f(t). For ABSDA day # betwen 2252 (May 1, 1984) and 4412 (Mar. 31, 1990):

```
N2O(t)(SIO-98)= N2O(t)(SIO-93)*1.0058/f(t)
```

where:
```
f(t) = a0 + a1*t + a2*t^(2) + a3*t^(3) + a4*t^(4)
t = (ABSDA-3227)/365.25
a0:  1.00664
a1: -0.56994e-3
a2: -0.65398e-3
a3:  0.13083e-3
a4: -0.20742e-4
```

ABSDA is the 3rd column in the original data file, 3227 corresponds to Jan. 1, 1987. Please be aware the f(t) function is based on day (e.g. ABSDA), all data within the same date is adjusted by the same f(t) (not varying with different hours, minutes).

## Repeatability
I haven't found a good single source of information on repeatabilties, so here are some notes from the early papers.

Prinn et al. (1983) estimates monthly variances (so repeatability + variability) of ~5% for CH3CCl3, ~2% for CFC-11, CFC-12, and CCl4, and 0.5% for N2O.

Cunnold et al. (1983) notes repeat measurement of standards with a repeatability of 0.5% for CFC-11, although 1% discontinuities for CFC-11 and CFC-12 found in Cunnold et al. (1986) due to tank changes.

Prinn et al. (2000) states *AGAGE* repeatabilities of 0.07% for CFC-11, 0.04% for CFC-12, 0.15% for CFC-113, 0.32% for CH3CCl3, 0.16% for CCl4 and 0.03% for N2O. States that this is improved over GAGE by 1 - 3 times at stations with good environmental control (e.g., CGO), and 2 - 6 times for other stations (e.g. RPB and SMO). A factor of 3 or 6 increase over AGAGE would give (for GAGE):
- CFC-11: 0.21% - 0.42%
- CFC-12: 0.12% - 0.24%
- CFC-113: 0.45% - 0.9%
- CH3CCl3: 1% - 2%
- CCl4: 0.5% - 1%

N2O repeatability is 0.35% for ALE and 0.13% for GAGE (Prinn et al., 1990).

CH4 GAGE repeatability is 0.3% at CGO to 0.6% at MHD (Cunnold et al., 2002).

Based on these bits of information, I've estimated repeatabilties as follows:

|Species |ALE repeatability (%)|GAGE repeatability (%)|Notes|
|--------|----|----|---|
|CFC-11  |1   |0.4 |GAGE value is at the upper end of the Prinn (2000) range inferred above. ALE is conservatively the 1% tank discontinuities mentioned in Cunnold (1986), and half the monthly variability in Prinn (1983).|
|CFC-12  |1   |0.2 |As for CFC-11.|
|CFC-113 |1   |1   |Used the Prinn (2000) inferred value for both.|
|CH3CCl3 |2   |2   |As for CFC-113.|
|CCl4    |1   |1   |As for CFC-113.|
|N2O     |0.35|0.13|Prinn (1990).|
|CH4     |-   |0.6 |Cunnold (2002) upper range applied to all.|

If you propose a change to these values, make sure they are consistent with ```ale_gage_species.json```.

## Scale uncertainty

SIO-93 scale uncertainty estimated as 0.8% for CFC-11 and 0.5% for CFC-12 (Cunnold et al., 1994)

## Site locations

Prinn et al. (2000) has Adrigole inlet 52m above sea level and Cape Mears, Oregon inlet 32m above sea level. Have assumed inlets were 10m high and adjusted ```inlet_base_height``` attribute accordingly.

## References
Cunnold, D. M., Prinn, R. G., Rasmussen, R. A., Simmonds, P. G., Alyea, F. N., Cardelino, C. A., Crawford, A. J., Fraser, P. J., and Rosen, R. D.: The Atmospheric Lifetime Experiment 3. Lifetime Methodology and Application to Three Years of CFCl3 Data, Journal of Geophysical Research, 88, 8379–8400, https://doi.org/10.1029/JC088iC13p08379, 1983.

Cunnold, D. M., Prinn, R. G., Rasmussen, R. A., Simmonds, P. G., Alyea, F. N., Cardelino, C. A., Crawford, A. J., Fraser, P. J., and Rosen, R. D.: Atmospheric lifetime and annual release estimates for CFCl 3 and CF 2 Cl 2 from 5 years of ALE data, J. Geophys. Res., 91, 10797, https://doi.org/10.1029/JD091iD10p10797, 1986.

Cunnold, D. M., Fraser, P. J., Weiss, R. F., Prinn, R. G., Simmonds, P. G., Miller, B. R., Alyea, F. N., and Crawford, A. J.: Global trends and annual releases of CCI3F and CCI2F2 estimated from ALE/GAGE and other measurements from July 1978 to June 1991, Journal of Geophysical Research, 99, 1107–1126, https://doi.org/10.1029/93JD02715, 1994.

Cunnold, D. M., Steele, L. P., Fraser, P. J., Simmonds, P. G., Prinn, R. G., Weiss, R. F., Porter, L. W., O’Doherty, S., Langenfelds, R. L., Krummel, P. B., Wang, H. J., Emmons, L., Tie, X. X., and Dlugokencky, E. J.: In situ measurements of atmospheric methane at GAGE/AGAGE sites during 1985–2000 and resulting source inferences, Journal of Geophysical Research, 107, 4225, https://doi.org/10.1029/2001JD001226, 2002.

Prinn, R. G., Simmonds, P. G., Rasmussen, R. A., Rosen, R. D., Alyea, F. N., Cardelino, C. A., Crawford, A. J., Cunnold, D. M., Fraser, P. J., and Lovelock, J. E.: The Atmospheric Lifetime Experiment 1. Introduction, Instrumentation, and Overview, Journal of Geophysical Research, 88, 8353–8367, https://doi.org/10.1029/JC088iC13p08353, 1983.

Prinn, R., Cunnold, D., Rasmussen, R., Simmonds, P., Alyea, F., Crawford, A., Fraser, P., and Rosen, R.: Atmospheric Emissions and Trends of Nitrous Oxide Deduced From 10 Years of ALE–GAGE Data, Journal of Geophysical Research, 95, 18369–18385, https://doi.org/10.1029/JD095iD11p18369, 1990.

Prinn, R. G., Weiss, R. F., Fraser, P. J., Simmonds, P. G., Cunnold, D. M., Alyea, F. N., O’Doherty, S., Salameh, P., Miller, B. R., Huang, J., Wang, R. H. J., Hartley, D. E., Harth, C., Steele, L. P., Sturrock, G., Midgley, P. M., and McCulloch,  a.: A history of chemically and radiatively important gases in air deduced from ALE/GAGE/AGAGE, Journal of Geophysical Research, 105, 17751–17792, https://doi.org/10.1029/2000JD900141, 2000.
