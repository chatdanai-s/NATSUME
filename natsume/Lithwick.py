# Functions that only exist in the Lithwick (2012) literature
import numpy as np
from .classes import ComplexEccentricities, TTVSineCurve
from .common import get_b, get_Db, get_alpha, get_NormalizedResonanceDistance, get_outerPeriods, get_innerPeriods

# Lithwick disturbing functions
def get_f(alpha, j: int):
    """
    Returns sum of Laplace coefficients with order-unity values (denoted f),
    as described in Table 3 of Lithwick et al. (2012).

    Args:
        alpha: Ratio of orbit semi-major axes, see natsume.common.get_alpha (float or array)
        j: Number denoting j:j-1 mean motion resonance of system (integer)

    Returns:
        f: As described in Table 3 of Lithwick et al. (2012). (float or array)
    """
    return -(j * get_b(alpha, j)) - (alpha/2 * get_Db(alpha, j, order=1))

def get_g(alpha, j: int):
    """
    Returns sum of Laplace coefficients with order-unity values (denoted g),
    as described in Table 3 of Lithwick et al. (2012).

    Args:
        alpha: Ratio of orbit semi-major axes, see natsume.common.get_alpha (float or array)
        j: Number denoting j:j-1 mean motion resonance of system (integer)

    Returns:
        g: As described in Table 3 of Lithwick et al. (2012). (float or array)
    """
    if j == 2:
        # -1/(2 * alpha^2) s.t. alpha ~ (1/2)**(2/3) for inner perturber,
        # -2 * alpha for outer perturber, which are identical.
        correction = -1.25992104989
    else:
        correction = 0
    return (j-0.5) * get_b(alpha, j-1) + (alpha/2 * get_Db(alpha, j-1, order=1)) + correction

# Weighted average of free eccentricities
def get_Zfree(f, g, z: ComplexEccentricities):
    """
    Returns weighted average of free eccentricities as described in equation 10 of Lithwick et al. (2012).

    Args:
        f: Sum of Laplace coefficients with order-unity values, see natsume.Lithwick.get_f (float or array)
        g: Sum of Laplace coefficients with order-unity values, see natsume.Lithwick.get_g (float or array)
        z: Complex eccentricities of two planets in the system (class ComplexEccentricies)

    Returns:
        Z_free: As described in equation 10 of Lithwick et al. (2012). (float or array)
    """
    return (f * z.inner_e * np.exp(1j * z.inner_periastron)) + \
           (g * z.outer_e * np.exp(1j * z.outer_periastron))


# Inversion functions
def LithwickOuterInversion(innerTTV: TTVSineCurve, innerPeriod: float,
                           j: int, z: ComplexEccentricities, outerPeriod='none'):
    """
    Returns mass of outer planet analytically computed from inner TTV sine curve,
    as described in equation 8 of Lithwick et al. (2012).

    Args:
        innerTTV: Inner TTV amplitude and TTV period, in days (class TTVSineCurve)
        innerPeriod: Orbital period of inner planet, in days (float)
        j: Number denoting j:j-1 mean motion resonance of system (integer)
        z: Complex eccentricities of two planets in the system (class ComplexEccentricies)
        outerPeriod: Period of outer planet, in days, if known ('none' or float)

    Returns:
        mu: Mass of outer planet with unit of host star mass (array or float).

        If outerPeriod == 'none', mu will be an array with two possible mass solutions.

        If outerPeriod is given, mu will be a float with one possible mass solution.
    """
    if outerPeriod == 'none':
        outerPeriods = get_outerPeriods(innerPeriod, innerTTV.superperiod, j, N=1)
    else:
        outerPeriods = outerPeriod

    alpha = get_alpha(innerPeriod, outerPeriods)
    f = get_f(alpha, j)
    g = get_g(alpha, j)
    Delta = get_NormalizedResonanceDistance(innerPeriod, outerPeriods, j, N=1)
    Zfree = get_Zfree(f, g, z)

    massRatio = np.pi * innerTTV.amplitude * j**(2/3) * (j-1)**(1/3) * np.abs(Delta) / innerPeriod / \
                np.abs(f + 1.5 * np.conj(Zfree) / Delta)
    return massRatio

def LithwickInnerInversion(outerTTV: TTVSineCurve, outerPeriod: float,
                           j: int, z: ComplexEccentricities, innerPeriod='none'):
    """
    Returns mass of inner planet analytically computed from outer TTV sine curve,
    as described in equation 9 of Lithwick et al. (2012).

    Args:
        outerTTV: Outer TTV amplitude and TTV period, in days (class TTVSineCurve)
        outerPeriod: Orbital period of outer planet, in days (float)
        j: Number denoting j:j-1 mean motion resonance of system (integer)
        z: Complex eccentricities of two planets in the system (class ComplexEccentricies)
        innerPeriod: Period of inner planet, in days, if known ('none' or float)

    Returns:
        mu: Mass of inner planet with unit of host star mass (array or float).

        If innerPeriod == 'none', mu will be an array with two possible mass solutions.

        If innerPeriod is given, mu will be a float with one possible mass solution.
    """
    if innerPeriod == 'none':
        innerPeriods = get_innerPeriods(outerPeriod, outerTTV.superperiod, j, N=1)
    else:
        innerPeriods = innerPeriod

    alpha = get_alpha(innerPeriods, outerPeriod)
    f = get_f(alpha, j)
    g = get_g(alpha, j)
    Delta = get_NormalizedResonanceDistance(innerPeriods, outerPeriod, j, N=1)
    Zfree = get_Zfree(f, g, z)

    massRatio = np.pi * outerTTV.amplitude * j * np.abs(Delta) / outerPeriod / \
                np.abs(-g + 1.5 * np.conj(Zfree) / Delta)
    return massRatio
