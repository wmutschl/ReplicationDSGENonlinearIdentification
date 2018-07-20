% Makes all model objects independent of names (parameters and variables
% are numbered) and computes analytical derivatives of f evaluated at steady state. 
%
% Inputs:
%       DSGE_orig: structure with original (symbolic) DSGE model
%       param: structure with all information about parameters
%       Settings: structure with Settings
%
% Calls:
%       anal_deriv_f: Computes Jacobian and Magnus-Neudecker Hessian of f
%
% Outputs:
%       DSGE_symbolic: structure with (symbolic) DSGE model independent of names
%       numbers: structure containing size of variables and parameters
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [DSGE_symbolic,numbers] = ModelPostCalculations(DSGE_orig,param,Settings)
fprintf('PUTTING MODEL INTO GENERAL FRAMEWORK...\n');
%% Precomputation
f=DSGE_orig.f; SelectMat = DSGE_orig.SelectMat; sig = DSGE_orig.sig; Sigma = DSGE_orig.Sigma; etatilde = DSGE_orig.etatilde; dfstudt = DSGE_orig.dfstudt;
states_0 = DSGE_orig.states_0; states_1 = DSGE_orig.states_1; 
shocks_1 = DSGE_orig.shocks_1; shocks_2 = DSGE_orig.shocks_2;
controls_1 = DSGE_orig.controls_1; controls_2 = DSGE_orig.controls_2;
numbers.nu = length(shocks_1); numbers.nx = length(states_0); numbers.nv = numbers.nu+numbers.nx;
numbers.ny = length(controls_1); numbers.nd = size(SelectMat,1);

%% Compute analytically Jacobian and Magnus-Neudecker-Hessian of f
fprintf('COMPUTING ANALYTICAL JACOBIAN AND MAGNUS-NEUDECKER-HESSIAN OF MODEL EQUATIONS...');
    [gra,hes] = anal_deriv_f(f,[states_1 shocks_2],controls_2,[states_0 shocks_1],controls_1,Settings.approx);
fprintf('FINISHED!\n');

%% Evaluate f, gra, and hes at steady-state, i.e. eliminate future variables
fprintf('EVALUATING MODEL EQUATIONS, GRADIENT, AND HESSIAN AT STEADY-STATE...');
    vars1 = transpose([states_0 shocks_1 controls_1]); 
    vars2 = transpose([states_1 shocks_2 controls_2]);
    f = subs(f, vars2,vars1);
    gra = subs(gra, vars2,vars1);
    hes = subs(hes, vars2,vars1); 
fprintf('FINISHED!\n');

%% Create vectors of numbered parameters and numbered variables
    DSGE_symbolic.vars=sym(zeros(length(vars1),1));
    for i=1:length(vars1)
        DSGE_symbolic.vars(i) = sym(['vars' num2str(i)]);
    end

%% Make expressions independent of names
fprintf('MAKE ALL EXPRESSIONS INDEPENDENT OF NAMES...');
    old = [param.names transpose(vars1)];
    new = [param.numbered transpose(DSGE_symbolic.vars)]; 
    DSGE_symbolic.f=subs(f,old,new); 
    DSGE_symbolic.gra=subs(gra,old,new); 
    DSGE_symbolic.hes=subs(hes,old,new); 
    DSGE_symbolic.Sigma=subs(Sigma,old,new);
    DSGE_symbolic.etatilde=subs(etatilde,old,new);
    DSGE_symbolic.sig = subs(sig,old,new);    
    DSGE_symbolic.SelectMat = subs(SelectMat,old,new);    
    DSGE_symbolic.dfstudt = subs(dfstudt,old,new);    
    DSGE_symbolic.distribution = DSGE_orig.distribution;
fprintf('FINISHED!\n');
