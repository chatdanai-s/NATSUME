# NATSUME
**N**ear-resonant **A**nalytic **T**TV **S**olver for **U**nknown **M**ass **E**stimates (NATSUME) for Python 3 (Work in progress!)

A python 3 module which aims to quickly estimate non-transiting exoplanet masses in possible near Mean Motion Resonance (MMR) solutions from approximately sinusoidal Transit Timing Variation (TTV) signals.

The TTV mass inversion estimations are based from Lithwick's model for 1st order near MMR (https://doi.org/10.1088/0004-637X/761/2/122) and Deck-Agol's model for higher order near MMRs (https://doi.org/10.3847/0004-637X/821/2/96).

Please note that the outputs will only serve as estimates to quickly gauge appropriate priors for more robust parameter fitting via n-body simulation codes, combined with methods such as nested sampling or MCMC. Moreover, the models assume a one-star two-planet system, low (but non-zero) orbital eccentricities and perfect coplanarity with respect to the line of sight.
