% Computes symbolic analytical derivatives of model objects w.r.t. param_identif and writes them to an m-file for numeric evaluation.
%
% Inputs:
%       DSGE_Model: structure containing all information about model
%       approx: order of approximation
%
% Calls: anal_deriv_param_print2f to print analytical derivatives of model objects into a script file
%
% Outputs:  symbolic derivatives of inputs written into a file
%
% Based upon Schmitt-Grohé & Uribe (2013): Implementing Iskrev's Identifiabilty Test (http://www.columbia.edu/~mu2166/iskrev/iskrev.html)
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function anal_deriv_param(DSGE_Model,approx)
model=DSGE_Model.shortname; 
f = DSGE_Model.symbolic.f; 
gra = DSGE_Model.symbolic.gra;
hes = DSGE_Model.symbolic.hes; 
Sigma = DSGE_Model.symbolic.Sigma;
etatilde = DSGE_Model.symbolic.etatilde;
sig = DSGE_Model.symbolic.sig; 
dfstudt = DSGE_Model.symbolic.dfstudt;
SS = DSGE_Model.symbolic.SS;
param_identif = DSGE_Model.param.identif; 
vars = DSGE_Model.symbolic.vars;


fprintf('COMPUTING SYMBOLIC ANALYTICAL DERIVATIVES W.R.T. PARAMETERS\n')
%% CONSTRUCT DSS_Dparam
    fprintf('    - Steady State...')
    DSS_Dparam = jacobian(SS,param_identif);
    fprintf('Finished!\n')
%% CONSTRUCT Dvariables_numbered_Dparam BY TERMS
    fprintf('    - Model equations...')
    Df_Donlyparam = jacobian(f,param_identif);  
    Df_Dvars = jacobian(f, vars);
    fprintf('Finished!\n')
%% CONSTRUCT Dgra_Dparam BY TERMS
    fprintf('    - Gradient...')
    Dgra_Donlyparam = jacobian(gra(:),param_identif);
    Dgra_Dvars = jacobian(gra(:),vars);
    fprintf('Finished!\n')
%% CONSTRUCT Dhes_Dparam BY TERMS    
    fprintf('    - Hessian...')
    Dhes_Donlyparam = jacobian(hes(:), param_identif);
    Dhes_Dvars = jacobian(hes(:), vars);
    fprintf('Finished!\n')
%% CONSTRUCT DSigma_Dparam, Detatilde_Dparam and Dsig_Dparam
    fprintf('    - Objects related to shocks and measurement errors...')    
    DSigma_Dparam = jacobian(Sigma(:),param_identif);    
    Dsig_Dparam = jacobian(sig,param_identif);
    Dsigetatilde_Dparam = jacobian(sig*etatilde(:),param_identif);
    Detatilde_Dparam = jacobian(etatilde(:),param_identif);
    if strcmp(DSGE_Model.symbolic.distribution,'Student-t')
        Ddfstudt_Dparam = jacobian(dfstudt,param_identif);
    else
        Ddfstudt_Dparam = [];
    end
    fprintf('Finished!\n')
%% PRINT to file
fprintf('SAVING SYMBOLIC ANALYTICAL DERIVATIVES INTO M-FILE FOR FURTHER EVALUATION...')
anal_deriv_param_print2f(model,DSGE_Model.spec,approx,DSS_Dparam,Df_Donlyparam,Df_Dvars,Dgra_Donlyparam,Dgra_Dvars,Dhes_Donlyparam,Dhes_Dvars,DSigma_Dparam,Dsig_Dparam,Dsigetatilde_Dparam,Detatilde_Dparam,Ddfstudt_Dparam);
fprintf('Finished!\n')
