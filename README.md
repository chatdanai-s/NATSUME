# NATSUME  
**N**ear-resonant **A**nalytic **T**TV **S**olver for **U**nknown **M**ass **E**stimates (NATSUME) for Python 3.

A python 3 module which aims to quickly estimate non-transiting exoplanet masses in possible near Mean Motion Resonance (MMR) scenarios from approximately sinusoidal Transit Timing Variation (TTV) signals.

TTV mass inversions are transcribed from Lithwick's model for 1st order near MMR (https://doi.org/10.1088/0004-637X/761/2/122) and Deck-Agol's model for higher order near MMRs (https://doi.org/10.3847/0004-637X/821/2/96).

Installation
=====
Install from pypi:
```
pip install natsume-ttv
```


Usage Overview
=====

`natsume` is to be used with the following steps:

1. Define the orbital elements (orbital periods, eccentricities, longitudes of periastron), TTV signal parameters (amplitudes, superperiods), and MMR scenario (e.g. `'2:1'`, `'5:3'`) of the planet pair. All units are in days or degrees.  
2. Encode complex eccentricity and sinusoidal TTV information into the `natsume.classes.ComplexEccentricities` and `natsume.classes.TTVSineCurve` objects using `natsume.get_ComplexEccentricities` and `natsume.get_TTVSineCurve` functions.
3. Estimate inner or outer exoplanet masses using `natsume.EstimateInnerMass` or `natsume.EstimateOuterMass` functions.

TTV amplitude and superperiod characterization is left to the user, as multiple valid approaches exist.

Usage Example
=====

To illustrate the steps above, suppose the following planetary parameters for a pair in the Kepler-32 system

```python
# Planetary parameters for Kepler-32b/c
Pb = 5.901    # Inner period (days; Lithwick+12, Table 1)
Pc = 8.752    # Outer period (days; Lithwick+12, Table 1)
Vb = 0.0062   # Inner TTV Amplitude (days; Lithwick+12, Table 1)
Vc = 0.0077	  # Outer TTV Amplitude (days; Lithwick+12, Table 1)
Pttvb = 267.1 # TTV "superperiod" for planet b (days; Holczer+16, Table 5)
Pttvc = Pttvb # TTV "superperiod" for planet c (Assume identical superperiod
              # as inner TTV as Holczer+16 did not report the period)             
mmr = '3:2'   # MMR Scenario

# Eccentricity and arguments of periastron (deg)
# Calculated from Fabrycky+12, Table 6 eccentricity vectors
eb = 0.003822
ec = 0.002384
wb = -166.56
wc = +43.98

Mstar = 0.49  # Host star mass (solar masses; Fabrycky+12, Table 6)
```

We then encode system complex eccentricity and sinusoidal TTV informations into the `ComplexEccentricities` object with `natsume.get_ComplexEccentricities` and `TTVSineCurve` object with `natsume.get_TTVSineCurve`

```python
import natsume

# Build sinusoidal TTV object for inner TTV
TTVb = natsume.get_TTVSineCurve(amplitude=Vb, superperiod=Pttvb)
TTVc = natsume.get_TTVSineCurve(amplitude=Vc, superperiod=Pttvc)

# Build complex eccentricity object for Kepler-32b/c pair
z = natsume.get_ComplexEccentricities(e1=eb, w1=wb, e2=ec, w2=wc)
```

Finally, we calculate exoplanet masses using analytic sinusoidal TTV models

```python
# Estimate outer planet mass relative to the host star
mu_c = natsume.EstimateOuterMass(
	innerTTV=TTVb,
	inner_period=Pb,
	mmr=mmr,
	eccentricity=z,
	outer_period=None
)
# Estimate inner planet mass relative to the host star
mu_b = natsume.EstimateInnerMass(
	outerTTV=TTVc,
	outer_period=Pc,
	mmr=mmr,
	eccentricity=z,
	inner_period=None
)
```

However, through above code, `mu_c` and `mu_b` will be `numpy` arrays containing two possible mass solutions. This is because the perturbing planet's period can be unknown if it is non-transiting, and two periods are possible given the definition of the TTV superperiod (see Lithwick's eqn. 5). If the perturbing planets' orbital periods are already known, we can skip the calculation of perturbing planet's orbital period entirely, and replace the `None` arguments with appropriate values

```python
# Estimate outer planet mass relative to the host star
mu_c = natsume.EstimateOuterMass(
	innerTTV=TTVb,
	inner_period=Pb,
	mmr=mmr,
	eccentricity=z,
	outer_period=Pc
)
# Estimate inner planet mass relative to the host star
mu_b = natsume.EstimateInnerMass(
	outerTTV=TTVc,
	outer_period=Pc,
	mmr=mmr,
	eccentricity=z,
	inner_period=Pb
)
```

To which `mu_c` and `mu_b` will be floats

Addendum
=====

Conversion from `mu_c` and `mu_b` which are relative to host stellar masses to `m_c` and `m_b`in Earth masses can also be done with `astropy.units`
```python
# Convert planet masses relative to host star to Earth masses
from astropy import units as u
m_b = (mu_b * Mstar*u.M_sun).to(u.M_earth).value
m_c = (mu_c * Mstar*u.M_sun).to(u.M_earth).value
```

Nominal mass calculations, where orbits are assumed with zero eccentricity, can be done by using the following variable `z0` in the `eccentricity` arguments of mass calculation functions
```python
# Build complex eccentricity object for zero eccentricity
z0 = natsume.get_ComplexEccentricities()
```
to which nominal mass calculations are only possible in near first-order resonance scenarios.
