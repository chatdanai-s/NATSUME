import numpy as np

from natsume.Lithwick import LithwickOuterInversion, LithwickInnerInversion
from natsume.DeckAgol import DeckAgolOuterInversion, DeckAgolInnerInversion
from natsume.common import MMR


### Classes
# Complex eccentricitries containing the eccentricity and
# longitude of periastron for inner and outer planets
class ComplexEccentricitries:
    def __init__(self,
                 inner_e=0, inner_periastron=0,
                 outer_e=0, outer_periastron=0):
        try:
            self.inner_e = float(inner_e)
            self.inner_periastron = float(inner_periastron)
            self.outer_e = float(outer_e)
            self.outer_periastron = float(outer_periastron)
        except (TypeError, ValueError):
            raise TypeError(f'ComplexEccentricitries arguments must be a float or integer.')

    @property
    def arr(self):
        return np.array([self.inner_e, self.inner_periastron,
                         self.outer_e, self.outer_periastron])

# TTV Sine curve with amplitude in minutes and superperiod in days
class TTVSineCurve:
    def __init__(self, amplitude, superperiod):
        try:
            self.amplitude = float(amplitude)
            self.superperiod = float(superperiod)
        except (TypeError, ValueError):
            raise TypeError(f'TTVSineCurve arguments must be a float or integer.')
        
    @property
    def arr(self):
        return np.array([self.amplitude, self.superperiod])


### Mass estimation functions
def EstimateOuterMass(innerTTV: TTVSineCurve, inner_period: float, mmr: str,
                      eccentricity=ComplexEccentricitries):
    j, N = MMR(mmr)

    if N == 1:
        mass = LithwickOuterInversion(innerTTV, inner_period, j, eccentricity)
    else:
        mass = DeckAgolOuterInversion(innerTTV, inner_period, j, N, eccentricity)

    return mass

def EstimateInnerMass(outerTTV: TTVSineCurve, outer_period: float, mmr: str,
                      eccentricity=ComplexEccentricitries):
    j, N = MMR(mmr)
    
    if N == 1:
        mass = LithwickInnerInversion(outerTTV, outer_period, j, eccentricity)
    else:
        mass = DeckAgolInnerInversion(outerTTV, outer_period, j, N, eccentricity)

    return mass
