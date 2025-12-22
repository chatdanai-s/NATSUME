# Near-resonant Analytic TTV Solver for Unknown Mass Estimation
from .classes import TTVSineCurve, ComplexEccentricities
from .common import get_MMR
from .Lithwick import LithwickOuterInversion, LithwickInnerInversion
from .DeckAgol import DeckAgolOuterInversion, DeckAgolInnerInversion

# Create classes (Cannot define classes here due to circular import)
def get_ComplexEccentricities(inner_e=0, inner_periastron=0,
                              outer_e=0, outer_periastron=0):
    """
    Returns ComplexEccentricities object for two planets in the system.

    Output used as input in functions requiring eccentricies and longitudes of periastron.

    Args:
        inner_e: Inner orbital eccentricity, must be at least 0 (float)
        inner_periastron: Inner longitude of periastron, in degrees (float)
        outer_e: Outer orbital eccentricity, must be at least 0 (float)
        outer_periastron: Outer longitude of periastron, in degrees (float)

    Returns:
        Complex eccentricities of two planets in the system (class ComplexEccentricities)
    """
    return ComplexEccentricities(inner_e, inner_periastron, outer_e, outer_periastron)

def get_TTVSineCurve(amplitude: float, superperiod: float):
    """
    Returns TTVSineCurve object of a planet in the system.
    Does not specify if inner TTV or outer TTV.

    Output used as input in functions requiring sinusoidal TTV curves.

    Args:
        amplitude: TTV amplitude in days, must be at least 0 (float)
        superperiod: TTV superperiod in days, must be at least 0 (float)

    Returns:
        TTV sine curve representation of a planet (class TTVSineCurve)
    """
    return TTVSineCurve(amplitude, superperiod)


### Mass estimation functions -- returns array of 2 if outer_period='none'
def EstimateOuterMass(innerTTV: TTVSineCurve, inner_period: float, mmr: str,
                      eccentricity=ComplexEccentricities(), outer_period='none'):
    """
    Returns mass of outer planet analytically computed from inner TTV sine curve,
    as described in equation 8 of Lithwick et al. (2012) or equation 18/42 of Deck and Agol (2016).

    Args:
        innerTTV: Inner TTV amplitude and TTV period, in days (class TTVSineCurve)
        inner_period: Orbital period of inner planet, in days (float)
        mmr: String denoting j:j-N mean motion resonance of system e.g. '2:1', '3:2', '3:1', etc. (string)
        eccentricity: Complex eccentricities of two planets in the system (class ComplexEccentricies)
        outer_period: Period of outer planet, in days, if known ('none' or float)

    Returns:
        mu: Mass of outer planet with unit of host star mass (array or float).

        If outer_period == 'none', mu will be an array with two possible mass solutions.

        If outer_period is given, mu will be a float with one possible mass solution.
    """
    j, N = get_MMR(mmr)
    if N == 1:
        mass = LithwickOuterInversion(innerTTV, inner_period, j,
                                      eccentricity, outerPeriod=outer_period)
    else:
        mass = DeckAgolOuterInversion(innerTTV, inner_period, j, N,
                                      eccentricity, outerPeriod=outer_period)
    return mass

def EstimateInnerMass(outerTTV: TTVSineCurve, outer_period: float, mmr: str,
                      eccentricity=ComplexEccentricities(), inner_period='none'):
    """
    Returns mass of inner planet analytically computed from outer TTV sine curve,
    as described in equation 9 of Lithwick et al. (2012) or equation 19/43 of Deck and Agol (2016).

    Args:
        outerTTV: Outer TTV amplitude and TTV period, in days (class TTVSineCurve)
        outer_period: Orbital period of outer planet, in days (float)
        mmr: String denoting j:j-N mean motion resonance of system e.g. '2:1', '3:2', '3:1', etc. (string)
        eccentricity: Complex eccentricities of two planets in the system (class ComplexEccentricies)
        inner_period: Period of inner planet, in days, if known ('none' or float)

    Returns:
        mu: Mass of inner planet with unit of host star mass (array or float).

        If inner_period == 'none', mu will be an array with two possible mass solutions.

        If inner_period is given, mu will be a float with one possible mass solution.
    """
    j, N = get_MMR(mmr)
    if N == 1:
        mass = LithwickInnerInversion(outerTTV, outer_period, j,
                                      eccentricity, innerPeriod=inner_period)
    else:
        mass = DeckAgolInnerInversion(outerTTV, outer_period, j, N,
                                      eccentricity, innerPeriod=inner_period)
    return mass
