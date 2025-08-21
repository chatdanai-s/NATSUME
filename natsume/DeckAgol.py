import numpy as np

from natsume.classes import ComplexEccentricities, TTVSineCurve

# Functions that only exist in the Deck-Agol (2016) literature
# Inversion functions
def DeckAgolOuterInversion(innerTTV: TTVSineCurve, inner_period: float,
                           j: int, N: int, eccentricity: ComplexEccentricities):
    pass

def DeckAgolInnerInversion(outerTTV: TTVSineCurve, outer_period: float,
                           j: int, N: int, eccentricity: ComplexEccentricities):
    pass