# Functions that only exist in the Lithwick (2012) literature
import numpy as np
from natsume import ComplexEccentricities, TTVSineCurve
from natsume.common import b, db_da, alpha, NormalizedResonanceDistance, OuterPeriods

def f(alpha, j: int):
    return -(j * b(alpha, j)) - (alpha/2 * db_da(alpha, j))

def g(alpha, j: int):
    return (j-0.5) * b(alpha, j-1) + (alpha/2 * db_da(alpha, j))

# Weighted average of free eccentricities
def Zfree(f, g, z: ComplexEccentricities):
    return (f * z.inner_e * np.exp(0+1j * z.inner_periastron)) + \
           (g * z.outer_e * np.exp(0+1j * z.outer_periastron))


# Inversion functions
def LithwickOuterInversion(innerTTV: TTVSineCurve, innerPeriod: float,
                           j: int, z: ComplexEccentricities):
    outerPeriods = OuterPeriods(innerPeriod, TTVSineCurve.superperiod, j, N=1)
    alpha = alpha(innerPeriod, outerPeriods)
    f = get_f()
    g = get_g()
    Delta = NormalizedResonanceDistance(innerPeriod, outerPeriods, j, N=1)
    Zfree = Zfree(f, g, z)
    massRatio = np.pi * TTVSineCurve.amplitude * j**(2/3) * (j-1)**(1/3) * Delta / innerPeriod / \
                np.abs(f + 1.5 * np.conj(Zfree) / Delta)
    return massRatio

def LithwickInnerInversion(outerTTV: TTVSineCurve, outerPeriod: float,
                           j: int, z: ComplexEccentricities):
    pass