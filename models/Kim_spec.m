% Model by Kim (2003) - Functional equivalence between intertemporal and
% multisectoral investment adjustment costs, in: Journal of Economic
% Dynamics and Control 27.4, pp. 533-549
% Output: 
%       spectext: text with specification of model
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

function spectext = Kim_spec
%% Model specification (spec)
%    0: no measurement errors, shock on a is iid normal
%    1: measurement error on all observables (normal), shock on a is normal
%    2: no measurement errors, shock on a is t-distributed
spectext = [...
    '0: no measurement errors, shock on a is normal\n '...
    '1: measurement error on all observables (normal), shock on a is normal\n'...
    '2: no measurement errors, shock on a is t-distributed'];
end