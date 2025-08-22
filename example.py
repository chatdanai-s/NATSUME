import natsume
from astropy import units as u

# # Kepler-32, all time unit in days (Lithwick et al. 2012)
# # Expected solution: 7.59 Earths for planet c
# Pb = 5.901
# Pc = 8.752 
# Vb = 0.0062
# PTTV = 1/abs(3/Pc - 2/Pb)
# mmr = '3:2'
# Mstar = 0.49
# z = natsume.get_ComplexEccentricities() # Zero eccentricity

# Kepler-46, all time unit in days (Saad-Olivera et al. 2017)
# Expected solution: 115 Earths for planet c
Pb = 33.648
Pc = 57.325
Vb = (70 * u.min).to(u.day).value
PTTV = 1/abs(5/Pc - 3/Pb)
mmr = '5:3'
Mstar = 0.902
z = natsume.get_ComplexEccentricities(0.0321, 264.2, 0.0354, 294.16) # Zero eccentricity

TTVb = natsume.get_TTVSineCurve(amplitude=Vb, superperiod=PTTV)

mu = natsume.EstimateOuterMass(innerTTV=TTVb, inner_period=Pb, mmr=mmr, eccentricity=z, outer_period=Pc)
m = (mu * Mstar*u.M_sun).to(u.M_earth).value

print(f'Estimated outer planet mass: {m} Earths')
