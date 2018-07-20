% Calculates steady-state either analytically or numerically for the 
% Model by Kim (2003) - Functional equivalence between intertemporal and
% multisectoral investment adjustment costs, in: Journal of Economic
% Dynamics and Control 27.4, pp. 533-549
% 
% Inputs: 
%       param: either symbolic or numerical values of all parameters, if symbolic independent of names, i.e. numbered
%       param_names: names of all parameters
%       spec: specification of model
% Outputs: 
%       SS: Numerical or symbolic expression for steady-state of the DSGE model
%
% Modified February 20, 2014 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

function [SS] = Kim_SS(param,param_names,spec)
%% Model specification (spec)
%    0: no measurement errors, shock on a is iid normal
%    1: measurement error on all observables (normal), shock on a is normal
%    2: no measurement errors, shock on a is t-distributed

% Get and evaluate parameter names and/or values
for i=1:length(param)   
    eval([char(param_names(i)) '= param(' num2str(i) ');'])
end

%% Declare steady state
% Auxiliary parameter
s=betae*delta*alph/(1-betae+delta*betae);

%% Define Steady-state of the model
% States
a=1;
k=(delta/s/a)^(1/(alph-1));
states = [k a];
% Shocks
e_a = 0; shocks = e_a;
% measurement errors
v_c=0; v_i=0;
if spec == 0
    meas_err = [];
elseif spec == 1
    meas_err = [v_c v_i];
elseif spec == 2
    meas_err = [];
end
 
% Controls
i=delta*k; 
c=(((a*k^alph)^(1+thet)-s*(i/s)^(1+thet))/(1-s))^(1/(1+thet))*(1-s);
controls = [c i];
% Observables
obs_c=c; obs_i=i;
observables = [obs_c obs_i];

% Define steady-state % DO NOT CHANGE!
SS = transpose([states shocks meas_err controls observables]);

