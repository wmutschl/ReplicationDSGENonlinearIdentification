% Parameter names and values for the 
% Model by Kim (2003) - Functional equivalence between intertemporal and
% multisectoral investment adjustment costs, in: Journal of Economic
% Dynamics and Control 27.4, pp. 533-549
%
%
% Inputs:
%       spec: Specification of model
% Outputs: 
%       params: structure containing information about parameters, values, priors etc.
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

function params = Kim_param(spec)

%% Model specification (spec)
%    0: no measurement errors, shock on a is iid normal
%    1: measurement error on all observables (normal), shock on a is normal
%    2: no measurement errors, shock on a is t-distributed

%% Set parameter names, keep same order as in model.m
% Prior types, please use UPPER Cases: 
% NORMAL, par1 is mean, par2 is std. deviation
% BETA, par1 is alpha, par2 ist beta
% GAMMA, par1 is shape parameter, par2 is scale parameter
% INVGAMMA, par1 is related to mean (location), par2 to tail to the right see Zellner (1971)
% UNIFORM, par1 is lower bound, par2 is upper bound
%

% Make parameters symbolic
syms alph betae delta as thet phi rho_a sigma_a sigma_c sigma_i dumpy dfstudt;


% each line has following structure
%         symbolic_name, latex_name, default_value, prior_type, prior_par1, prior_par2, bound_type, bound_low, bound_up, bound_parameter_c
params{1} = {alph,      '\alpha',      '0.60',       'GAMMA',      '0.60',     '0.30',      'R^+',     '1e-5',    '1',     '1'};
params{2} = {betae,     '\beta',       '0.99',       'UNIFORM',    '0.95',     '0.9999',    '[a,b)',   '0.9',    '0.99999','1'};
params{3} = {delta,     '\delta',      '0.0125',     'UNIFORM',    '0.01',     '0.02',      '[a,b)',   '0.01',    '0.02',   '1'};
params{4} = {thet,      '\theta',      '1',          'NORMAL',     '1.00',     '0.50',      'R',       '-5',      '5',      '1'};
params{5} = {phi,       '\phi',        '2',          'NORMAL',     '2.00',     '0.50',      'R',       '-5',      '5',      '1'};
params{6} = {rho_a,     '\rho_a',      '0.7',        'BETA',       '0.50',     '0.20',      '[a,b)',   '1e-5',    '0.99999','1'};
params{7} = {sigma_a,   '\sigma_a',    '0.5',        'INVGAMMA',   '0.50',     '4.00',      'R^+',       '1e-8',    '5',      '1'};
params{8} = {dumpy,     'dumpy'        '1',          'NORMAL',     '1.00',     '0.2',        'R',      '-5',      '5',       '1'};
meas_c = {sigma_c,      '\sigma_c',    '0.5',          'INVGAMMA',   '0.50',     '4.00',      'R^+',     '1e-5',    '10',     '1'};
meas_i = {sigma_i,      '\sigma_i',    '0.5',          'INVGAMMA',   '0.50',     '4.00',      'R^+',     '1e-5',    '10',     '1'};
df =     {dfstudt,            'df_t',  '10',         'UNIFORM',     '8',        '20',        '[a,b)',    '8',    '30',      '1'};
if spec == 1
    params{9} = meas_c;
    params{10} = meas_i;
elseif spec == 2
    params{9} = df;
end


