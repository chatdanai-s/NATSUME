# Functions that only exist in the Deck-Agol (2016) literature
import numpy as np

from natsume.classes import ComplexEccentricities, TTVSineCurve
from natsume.common import get_outerPeriods, get_innerPeriods, get_alpha, get_NormalizedResonanceDistance
from natsume.common import get_b, get_db_da, get_d2b_da2

# Laplace coefficients (2nd order)
def get_g_j45(alpha, j):
    return 1/8 * ((-5*j + 4*j**2) * get_b(alpha, j) + \
                  (4*j - 2) * alpha * get_db_da(alpha, j) + \
                  alpha**2 * get_d2b_da2(alpha, j))

def get_g_j49(alpha, j):
    return 1/4 * ((-2 + 6*j + 4*j**2) * get_b(alpha, j-1) - \
                  (4*j - 2) * alpha * get_db_da(alpha, j-1) - \
                  alpha**2 * get_d2b_da2(alpha, j-1))

def get_g_j53(alpha, j):
    if j == 3: # 27/8 * alpha where alpha ~ (1/3)**(2/3)
        correction = -1.6225307666
    else:
        correction = 0
    return 1/8 * ((2 - 7*j + 4*j**2) * get_b(alpha, j-2) + \
                  (4*j - 2) * alpha * get_db_da(alpha, j-2) + \
                  alpha**2 * get_d2b_da2(alpha, j-2)) + correction

# 2nd order AB functions -- Laplace coeffs will be separate args to reduce integration
def get_A1(g45, g49, g53, e1, e2, w1, w2):
    return (g45 * e1**2 * np.cos(2*w1)) + \
           (g53 * e2**2 * np.cos(2*w2)) + \
           (g49 * e1*e2 * np.cos(w1+w2))

def get_A2(g45, g49, g53, e1, e2, w1, w2):
    return -(g45 * e1**2 * np.sin(2*w1)) - \
            (g53 * e2**2 * np.sin(2*w2)) - \
            (g49 * e1*e2 * np.sin(w1+w2))

def get_B11(g45, g49, e1, e2, w1, w2):
    return (2 * g45 * e1 * np.cos(w1)) + (g49 * e2 * np.cos(w2))

def get_B12(g45, g49, e1, e2, w1, w2):
    return -(2 * g45 * e1 * np.sin(w1)) - (g49 * e2 * np.sin(w2))

def get_B21(g49, g53, e1, e2, w1, w2):
    return (2 * g53 * e2 * np.cos(w2)) + (g49 * e1 * np.cos(w1))

def get_B22(g49, g53, e1, e2, w1, w2):
    return -(2 * g53 * e2 * np.sin(w2)) - (g49 * e1 * np.sin(w1))


# Inversion functions
def DeckAgolOuterInversion(innerTTV: TTVSineCurve, innerPeriod: float,
                           j: int, N: int, eccentricity: ComplexEccentricities, outerPeriod='none'):
    if outerPeriod == 'none':
        outerPeriods = get_outerPeriods(innerPeriod, innerTTV.superperiod, j, N)
    else:
        outerPeriods = outerPeriod
    
    alpha = get_alpha(innerPeriod, outerPeriods)
    Delta = get_NormalizedResonanceDistance(innerPeriod, outerPeriods, j, N)
    e1, w1, e2, w2 = eccentricity.arr
    w1, w2 = np.deg2rad(w1), np.deg2rad(w2) # Uses rad in np.cos() and np.sin()

    if (e1 == 0) and (e2 == 0):
        raise ValueError('The Deck-Agol model does not provide physical zero-eccentricity mass solutions at N > 1.')

    if N == 2: # Case 2nd Order
        g45 = get_g_j45(alpha, j)
        g49 = get_g_j49(alpha, j)
        g53 = get_g_j53(alpha, j)

        A1 = get_A1(g45, g49, g53, e1, e2, w1, w2)
        A2 = get_A2(g45, g49, g53, e1, e2, w1, w2)
        B11 = get_B11(g45, g49, e1, e2, w1, w2)
        B12 = get_B12(g45, g49, e1, e2, w1, w2)

        massRatio = np.pi * innerTTV.amplitude * j**(2/3) * (j-N)**(1/3) * np.abs(Delta) / innerPeriod / \
                    np.sqrt((1.5 * A1 / Delta + B11)**2 + (1.5 * A2 / Delta + B12)**2)
    else:
        raise NotImplementedError('NATSUME does not yet support N = 3 and beyond near MMR.')
    
    return massRatio

def DeckAgolInnerInversion(outerTTV: TTVSineCurve, outerPeriod: float,
                           j: int, N: int, eccentricity: ComplexEccentricities, innerPeriod='none'):
    if innerPeriod == 'none':
        innerPeriods = get_innerPeriods(outerPeriod, outerTTV.superperiod, j, N)
    else:
        innerPeriods = innerPeriod

    alpha = get_alpha(innerPeriods, outerPeriod)
    Delta = get_NormalizedResonanceDistance(innerPeriods, outerPeriod, j, N=1)
    e1, w1, e2, w2 = eccentricity.arr
    w1, w2 = np.deg2rad(w1), np.deg2rad(w2) # Uses rad in np.cos() and np.sin()
    
    if (e1 == 0) and (e2 == 0):
        raise ValueError('The Deck-Agol model does not provide physical zero-eccentricity mass solutions at N > 1.')
    
    if N == 2: # Case 2nd order
        g45 = get_g_j45(alpha, j)
        g49 = get_g_j49(alpha, j)
        g53 = get_g_j53(alpha, j)

        A1 = get_A1(g45, g49, g53, e1, e2, w1, w2)
        A2 = get_A2(g45, g49, g53, e1, e2, w1, w2)
        B21 = get_B21(g49, g53, e1, e2, w1, w2)
        B22 = get_B22(g49, g53, e1, e2, w1, w2)

        massRatio = np.pi * outerTTV.amplitude * j * Delta / outerPeriod / \
                    np.sqrt((-1.5 * A1 / Delta + B21)**2 + (-1.5 * A2 / Delta + B22)**2)
    else:
        raise NotImplementedError('NATSUME does not yet support N = 3 and beyond near MMR.')
    
    return massRatio
