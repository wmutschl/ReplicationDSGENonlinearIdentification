% Parameter names and values for the 
% Model by An and Schorfheide (2007) - Bayesian Analysis of DSGE models, in: Econometric Reviews 26(2-4):113-172.
%
%
%
% Inputs:
%       spec: Specification of model
% Outputs: 
%       params: structure containing information about parameters, values, priors etc.
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)


function params = AnSchorfheide_param(spec)

%% Model specification (spec)
% 0: no measurement errors
% 1: measurement error on all observables
% 2: no measurement errors, all shocks are t-distributed

%% Set parameter names, keep same order as in model.m
% Prior types, please use UPPER Cases: 
% NORMAL, par1 is mean, par2 is std. deviation
% BETA, par1 is alpha, par2 ist beta
% GAMMA, par1 is shape parameter, par2 is scale parameter
% INVGAMMA, par1 is related to mean (location), par2 to tail to the right see Zellner (1971)
% UNIFORM, par1 is lower bound, par2 is upper bound
%

% Make parameters symbolic
syms tau phi psi1 psi2 rhoR rhog rhoz rA pA gamQ sigR sigg sigz nu cyst sig_YGR sig_INFL sig_INT dumpy dfstudt;

% each line has following structure
%         symbolic_name, latex_name, default_value, prior_type, prior_par1, prior_par2, bound_type, bound_low, bound_up, bound_parameter_c
params{1} = {tau,      '\tau',          '2.00',     'GAMMA',       '2.00',     '0.50',      'R^+',       '1e-5',    '10',     '1'};
params{2} = {phi,      '\phi',          '50',       'GAMMA',       '50',       '20',        'R^+',       '1e-5',    '100',    '1'};
params{3} = {psi1,     '\psi_1',        '1.50',     'GAMMA',       '1.50',     '0.25',      'R^+',       '1e-5',    '10',     '1'};
params{4} = {psi2,     '\psi_2',        '0.125',    'GAMMA',       '0.50',     '0.25',      'R^+',       '1e-5',    '10',     '1'};
params{5} = {rhoR,     '\rho_R',        '0.75',     'BETA',        '0.50',     '0.20',      '[a,b)',     '1e-5',    '0.99999','1'};
params{6} = {rhog,     '\rho_g',        '0.95',     'BETA',        '0.80',     '0.10',      '[a,b)',     '1e-5',    '0.99999','1'};
params{7} = {rhoz,     '\rho_z',        '0.90',     'BETA',        '0.66',     '0.15',      '[a,b)',     '1e-5',    '0.99999','1'};
params{8} = {rA,       'r^{(A)}',       '1.00',     'GAMMA',       '0.80',     '0.50',      'R^+',       '1e-5',    '10',     '1'};
params{9} = {pA,       'p^{(A)}',       '3.20',     'GAMMA',       '4.00',     '2.00',      'R^+',       '1e-5',    '20',     '1'};
params{10} = {gamQ,    '\gamma_Q',      '0.55',     'NORMAL',      '0.40',     '0.20',      'R',         '-5',      '5',      '1'};
params{11} = {sigR,    '\sigma_R',      '0.002',    'INVGAMMA',    '0.30',     '4.00',      'R^+',       '1e-8',    '5',      '1'};
params{12} = {sigg,    '\sigma_g',      '0.006',    'INVGAMMA',    '0.40',     '4.00',      'R^+',       '1e-8',    '5',      '1'};
params{13} = {sigz,    '\sigma_z',      '0.003',    'INVGAMMA',    '0.40',     '4.00',      'R^+',       '1e-8',    '5',      '1'};
params{14} = {nu,       '\nu',           '0.10',     'BETA',        '0.10',     '0.05',      '[a,b)',     '1e-5',    '0.99999','1'};
params{15} = {cyst,    'c/y',           '0.85',     'BETA',        '0.85',     '0.10',      '[a,b)',     '1e-5',    '0.99999','1'};
params{16} = {dumpy,     'dumpy'        '1',          'NORMAL',     '1.00',     '0.2',        'R',      '-5',      '5',       '1'};
if spec == 1
params{17} = {sig_YGR, '\sigma_{YGR}',  '0.23',     'INVGAMMA',    '0.30',     '4.00',      'R^+',       '1e-8',    '5',      '1'};
params{18} = {sig_INFL,'\sigma_{INFL}', '0.56',     'INVGAMMA',    '0.50',     '4.00',      'R^+',       '1e-8',    '5',      '1'};
params{19} = {sig_INT, '\sigma_{INT}',  '0.66',     'INVGAMMA',    '0.60',     '4.00',      'R^+',       '1e-8',    '5',      '1'};
end
if spec == 2
    params{17} =     {dfstudt,            'df_{t}',           '10',         'UNIFORM',     '8',        '20',        '[a,b)',    '8',    '35',      '1'};
end

