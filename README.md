# Replication files for Mutschler - Identification of DSGE models - The effect of second order approximation and pruning, Journal of Economic Dynamics and Control, Volume 56, July 2015, Pages 34-54, [DOI](http://dx.doi.org/10.1016/j.jedc.2015.04.007)

This is the documentation for the code and also contains some additional material

The paper establishes rank criteria for local identification given the pruned state-space representation in the fashion of Iskrev (2010) and Qu and Tkachenko (2012), also including higher-order moments, cumulants and polyspectra. It is shown that this may improve overall identification of a DSGE model via imposing additional restrictions on the moments and spectra.

In the Matlab code the user can choose in a graphical-user-interface between the models, the tests, which parameters to identify at which local point, analytical or numerical derivatives, and the order of approximation.

Since all procedures are model independent, other models can be easily included and tested as long as they can be represented in the same framework.

How to run:
- You will need Matlab's symbolic toolbox
- Make sure to be in the main directory
- Just run `identification_run.m` all options are asked via a GUI

The folder `nonidentification_curves` contains additional codes to compute nonidentification curves in the fashion of Qu and Tkachenko (2012) for robustness checks.

The folder `additional_material` contains additional material as mentioned in the paper.
