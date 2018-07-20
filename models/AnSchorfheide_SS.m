% Calculates steady-state either analytically or numerically for the 
% Model by An and Schorfheide (2007) - Bayesian Analysis of DSGE models, in: Econometric Reviews 26(2-4):113-172.
% 
% Inputs: 
%       param: either symbolic or numerical values of all parameters, if symbolic independent of names, i.e. numbered
%       param_names: names of all parameters
%       spec: specification of model
% Outputs: 
%       SS: Numerical or symbolic expression for steady-state of the DSGE model
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

function [SS] = AnSchorfheide_SS(param,param_names,spec)
%% Model specification (spec)
% 0: no measurement errors
% 1: measurement error on all observables
% 2: no measurement errors, all shocks are t-distributed

%% Do not change this block
% Get and evaluate parameter names and/or values
for i=1:length(param)   
    eval([char(param_names(i)) '= param(' num2str(i) ');'])
end

%% Define Steady-state of the model
% States
y = 0; R = 0; g = 0; z = 0; 
states = [y R g z];
% Shocks 
e_R = 0; e_g = 0; e_z = 0; 
% measurement errors
v_YGR=0; v_INFL=0; v_INT=0; 
if spec == 1
    meas_err = [v_YGR v_INFL v_INT];
else
    meas_err =[];
end
shocks = [e_R e_g e_z meas_err]; 
% Controls
c = 0; dy = 0; p = 0; 
YGR = gamQ; INFL = pA; INT = pA + rA + 4*gamQ; 
controls = [c dy p YGR INFL INT];

% Define steady-state % DO NOT CHANGE!
SS = transpose([states shocks controls]);

