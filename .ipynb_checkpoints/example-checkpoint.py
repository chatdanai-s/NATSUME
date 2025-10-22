import natsume
from astropy import units as u

# # Kepler-46, all time unit in days (Saad-Olivera et al. 2017)
# # Expected solution: 115 Earths for planet c
# Pb = 33.648
# Pc = 57.325
# Vb = (70 * u.min).to(u.day).value
# PTTV = 1/abs(5/Pc - 3/Pb)
# mmr = '5:3'
# Mstar = 0.902
# z = natsume.get_ComplexEccentricities(0.0321, 264.2, 0.0354, 294.16)

# # Kepler-32, all time unit in days (Taken from Lithwick et al. 2012)
# # Expected nominal solution: 7.59 pm 0.20 Earths for planet c
# # Actual solution: 6.0 pm 1.9 Earths for planet 
# Pb = 5.901
# Pc = 8.752
# Vb = 0.0062
# PTTV = 1/abs(3/Pc - 2/Pb)
# mmr = '3:2'
# Mstar = 0.49
# z = natsume.get_ComplexEccentricities() # Zero eccentricity

# Kepler-18-like/Figure 3, all time unit in days (Taken from Lithwick et al. 2012)
# Expected solution: 5.846 Earths for planet c
# Actual solution: ??? Earths for planet c
Pb = 7.642
Pc = 14.86
Vb = 0.002306
PTTV = 268.751482
mmr = '2:1'
Mstar = 0.972
z = natsume.get_ComplexEccentricities(0.05, 0, 0.05, 0)

TTVb = natsume.get_TTVSineCurve(amplitude=Vb, superperiod=PTTV)
mu = natsume.EstimateOuterMass(innerTTV=TTVb, inner_period=Pb, mmr=mmr, eccentricity=z)
m = (mu * Mstar*u.M_sun).to(u.M_earth).value

print(f'Estimated outer planet mass: {m} Earths')
