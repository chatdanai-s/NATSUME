import natsume
from astropy import units as u

# Kepler-32b, all time unit in days (Lithwick et al. 2012)
# Expected solution: 7.59 Earths
Pb = 5.901
Pc = 8.752 
Vb = 0.0062
PTTV = 1/abs(3/Pc - 2/Pb)
mmr = '3:2'

Mstar = 0.49

z = natsume.get_ComplexEccentricities() # Zero eccentricity
TTVb = natsume.get_TTVSineCurve(amplitude=Vb, superperiod=PTTV)

mu = natsume.EstimateOuterMass(innerTTV=TTVb, inner_period=Pb, mmr=mmr, eccentricity=z, outer_period=Pc)
m = (mu * Mstar*u.M_sun).to(u.M_earth).value

print(f'Estimated outer planet mass: {m} Earths')
