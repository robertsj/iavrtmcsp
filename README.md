# iavrtmcsp
Report and code for 22.106 final project Investigation of Automated Variance Reduction Techniques for Monte Carlo Shielding Problems

# Abstract

The Monte Carlo method is widely believed to be the most accurate method for solving problems in radiation transport.  Unfortunately, due to its very nature---following individual particle histories---certain classes of problems are particularly challenging for the method.  One such class of problems consist of so-called deep penetration shielding problems.  Because the purpose of a shield is to attenuate a particle population by several orders of magnitude, to use the Monte Carlo method requires a sufficient number of histories to ensure that the population, once attenuated, can still provide adequate statistics.  For deep penetration problems, the level of attenuation makes Monte Carlo prohibitively expensive.

To circumvent this issue, several approaches for *variance reduction* have been developed over the years.  Variance reduction techniques aim to modify (i.e., bias) in some manner the underlying physics in such a way that an *unbiased* solution with *lower* statistical error is found than an unbiased simulation using the same computational resources.  Haghighat and Wagner classify variance reduction techniques in three ways: *modified sampling methods* (e.g., source biasing, implicit capture), *population control methods* (e.g., geometry splitting/roulette, weight windows), and *semi-analytical methods* (e.g., point detectors and DXTRAN).  In this paper, we only analyze methods falling in the first two categories, namely (automated) approaches for source biasing, geometry splitting/roulette, and weight windows.

# Comments

This project is old (spring 2010) and was my first attempt to understand variance reduction for Monte Carlo transport based on deterministic methods.  I did not pursue further work in this area, but others at ORNL and elsewhere have continued the search for practical methods.  Two areas that might warrant more thought include adaptation of the largely source/detector framework to eigenvalue problems (in which we still care about the response of discrete sensors despite not having a fixed source) and time-dependent biasing for transient simulations.
