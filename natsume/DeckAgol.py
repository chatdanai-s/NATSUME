# Functions that only exist in the Deck-Agol (2016) literature
import numpy as np

from natsume.classes import ComplexEccentricities, TTVSineCurve
from natsume.common import get_outerPeriods, get_innerPeriods, get_alpha, get_NormalizedResonanceDistance
from natsume.common import get_b, get_Db

# Source: Dermott & Murray 1999, Appendix B
# Disturbing function coefficients (2nd order)
def get_g_j45(alpha, j): # k=N=2
    return 1/8 * ((-5*j + 4*j**2) * get_b(alpha, j) + \
                  (4*j - 2) * alpha * get_Db(alpha, j, order=1) + \
                  alpha**2 * get_Db(alpha, j, order=2))

def get_g_j49(alpha, j): # k=1
    return 1/4 * ((-2 + 6*j + 4*j**2) * get_b(alpha, j-1) - \
                  (4*j - 2) * alpha * get_Db(alpha, j-1, order=1) - \
                  alpha**2 * get_Db(alpha, j-1, order=2))

def get_g_j53(alpha, j): # k=0
    if j == 3:
        correction = -1.6225307666
        # -3/(8 * alpha^2) s.t. alpha ~ (1/3)**(2/3) for inner perturber,
        # -27/8 * alpha for outer perturber, which are identical!
    else:
        correction = 0
    return 1/8 * ((2 - 7*j + 4*j**2) * get_b(alpha, j-2) + \
                  (4*j - 2) * alpha * get_Db(alpha, j-2, order=1) + \
                  alpha**2 * get_Db(alpha, j-2, order=2)) + correction

# Disturbing function coefficients (3rd order)
def get_g_j82(alpha, j): # k=N=3
    return 1/48 * ((-26*j + 30*j**2 - 8*j**3) * get_b(alpha, j) + \
                   (-9 + 27*j - 12*j**2) * alpha * get_Db(alpha, j ,order=1) + \
                   (6 - 6*j) * alpha**2 * get_Db(alpha, j, order=2) - \
                   alpha**3 * get_Db(alpha,j , order=3))

def get_g_j83(alpha, j): # k=2
    return 1/16 * ((-9 + 31*j - 30*j**2 + 8*j**3) * get_b(alpha, j-1) + \
                   (9 - 25*j + 12*j**2) * alpha * get_Db(alpha, j-1, order=1) + \
                   (-5 + 6*j) * alpha**2 * get_Db(alpha, j-1, order=2) + \
                   alpha**3 * get_Db(alpha, j-1, order=3))

def get_g_j84(alpha, j): # k=1
    return 1/16 * ((8 - 32*j + 30*j**2 + 8*j**3) * get_b(alpha, j-2) + \
                   (-8 + 23*j - 12*j**2) * alpha * get_Db(alpha, j-2, order=1) + \
                   (4 - 6*j) * alpha**2 * get_Db(alpha, j-2, order=2) - \
                   alpha**3 * get_Db(alpha, j-2, order=3))

def get_g_j85(alpha, j): # k=0
    if j == 4:
        # -1/(3 * alpha^2) s.t. alpha ~ (1/4)**(2/3) for inner perturber,
        # -16/3 * alpha for outer perturber, which are identical!
        correction = -2.11653473596
    else:
        correction = 0
    return 1/48 * ((-6 + 29*j - 30*j**2 + 8*j**3) * get_b(alpha, j-3) + \
                   (6 - 21*j + 12*j**2) * alpha * get_Db(alpha, j-3, order=1) + \
                   (-3 + 6*j) * alpha**2 * get_Db(alpha, j-3, order=2) + \
                   alpha**3 * get_Db(alpha, j-3, order=3)) + correction

# Disturbing function coefficients (4th order)
def get_g_j90(alpha, j): # k=N=4
    return 1/384 * ((-206*j + 283*j**2 - 120*j**3 + 16*j**4) * get_b(alpha, j) + \
                    (-64 + 236*j - 168*j**2 + 32*j**3) * alpha * get_Db(alpha, j, order=1) + \
                    (48 - 78*j + 24*j**2) * alpha**2 * get_Db(alpha, j, order=2) + \
                    (-12 + 8*j) * alpha**3 * get_Db(alpha, j, order=3) + 
                    alpha**4 * get_Db(alpha, j, order=4))

def get_g_j91(alpha, j): # k=3
    return 1/96 * ((-64 + 238*j - 274*j**2 + 116*j**3 - 16*j**4) * get_b(alpha, j-1) + \
                   (64 - 206*j + 156*j**2 - 32*j**3) * alpha * get_Db(alpha, j-1, order=1) + \
                   (-36 + 69*j - 24*j**2) * alpha**2 * get_Db(alpha, j-1, order=2) + \
                   (10 - 8*j) * alpha**3 * get_Db(alpha, j-1, order=3) - \
                   alpha**4 * get_Db(alpha, j-1, order=4))

def get_g_j92(alpha, j): # k=2
    return 1/64 * ((52 - 224*j + 259*j**2 - 112*j**3 + 16*j**4) * get_b(alpha, j-2) + \
                   (-52 + 176*j - 144*j**2 + 32*j**3) * alpha * get_Db(alpha, j-2, order=1) + \
                   (26 - 60*j + 24*j**2) * alpha**2 * get_Db(alpha, j-2, order=2) + \
                   (-8 + 8*j) * alpha**3 * get_Db(alpha, j-2, order=3) + \
                   alpha**4 * get_Db(alpha, j-2, order=4))

def get_g_j93(alpha, j): # k=1
    return 1/96 * ((-36 + 186*j - 238*j**2 + 108*j**3 - 16*j**4) * get_b(alpha, j-3) + \
                   (36 - 146*j + 132*j**2 - 32*j**3) * alpha * get_Db(alpha, j-3, order=1) + \
                   (-18 + 51*j - 24*j**2) * alpha**2 * get_Db(alpha, j-3, order=2) + \
                   (6 - 8*j) * alpha**3 * get_Db(alpha, j-3, order=3) - \
                   alpha**4 * get_Db(alpha, j-3, order=4))

def get_g_j94(alpha, j): # k=0
    if j == 5:
        # -125/(384 * alpha^2) s.t. alpha ~ (1/5)**(2/3) for inner perturber,
        # -3125/384 * alpha for outer perturber, which are identical!
        correction = -2.78316397571
    else:
        correction = 0
    return 1/384 * ((24 - 146*j + 211*j**2 - 104*j**3 + 16*j**4) * get_b(alpha, j-4)
                    (-24 + 116*j - 120*j**2 + 32*j**3) * alpha * get_Db(alpha, j-4, order=1) + \
                    (12 - 42*j + 24*j**2) * alpha**2 * get_Db(alpha, j-4, order=2) + \
                    (-4 + 8*j) * alpha**3 * get_Db(alpha, j-4, order=3) + \
                    alpha**4 * get_Db(alpha, j-4, order=4)) + correction


# AB functions -- Laplace coeffs will be separate args to reduce integration iterations
# Will eventually be generalized to nth order
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
        raise NotImplementedError('NATSUME does not yet support TTVs of order N > 2 near MMR.')
    
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
        raise NotImplementedError('NATSUME does not yet support TTVs of order N > 2 near MMR.')
    
    return massRatio
