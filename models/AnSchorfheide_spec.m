% Model by An and Schorfheide (2007) - Bayesian Analysis of DSGE models, in: Econometric Reviews 26(2-4):113-172.
% Output: 
%       spectext: text with specification of model
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)
function spectext = AnSchorfheide_spec
%% Model specification (spec)
% 0: no measurement errors, all shocks are Gaussian
% 1: measurement error on all observables, all shocks are Gaussian
% 2: no measurement errors, all shocks are t-distributed
spectext = [...
    '0: no measurement errors \n '...
    '1: measurement error on all observables \n '...
    '2: no measurement errors, all shocks are t-distributed'];

end