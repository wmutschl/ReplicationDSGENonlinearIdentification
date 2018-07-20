% Evaluates numerically steady state, model objects and analytical Jacobians
%
% Inputs:
%       param_estim: vector with numerical values for all parameters
%       DSGE_Model: structure containing DSGE model independent of names
%       approx: order of approximation
%       Deriv_type: Analytical or Numerical derivatives
%
% Calls: 
%   -'modelshortname'_SS for numerical steady state evaluation
%   -'modelshortname'_num_eval for numerical evaluation of model objects
%   -'modelshortname'_anal_deriv for numerical evaluation of analytical derivatives of model objects
% 
% Outputs:
%       ngra: numerical evaluation of the gradient of f, evaluated at the steady-state
%       nf: numerical evaluation of the model equations f, evaluated at the steady-state
%       nhes: numerical evaluation of the Magnus-Neudecker-Hessian of f, evaluated at the steady-state
%       nSigma: numerical evaluation of the var/cov matrix of the shocks and measurement errors
%       netatilde: numerical evaluation of the matrix scaling the std. deviation of the shocks and measurement errors in the transition equation
%       nsig: numerical evaluation of the perturbation parameter
%       nSelectMat: numerical evaluation of the selection matrix
%       nSS: numerical evaluation of the steady state
%       ndfstudt: numerical evaluated degrees of freedom for student's t-distribution
%       nd,nx,ny,nu: Number of observables, states, control, shocks&measurement erros, nv=nx+nu
%       nDgra_Dparam: numerical evaluation of the gradient of f differentiated w.r.t param_identif, evaluated at the steady-state
%       nDhes_Dparam: numerical evaluation of the Magnus-Neudecker-Hessian of f differentiated w.r.t param_identif, evaluated at the steady-state
%       nDSigma_Dparam: numerical evaluation of the var/cov matrix of the shocks and measurement errors in the transition equation including the perturbation parameter differentiated w.r.t param_identif
%       nDetatilde_Dparam: numerical evaluation of the matrix scaling the std. deviation of the shocks and measurement errors in the transition equation without the perturbation parameter differentiated w.r.t param_identif
%       nDsig_Dparam: numerical evaluation of the perturbation parameter differentiated w.r.t param_identif
%       nDsigetatilde_Dparam: numerical evaluation of the derivative of sig*etatilde
%       nDSS_Dparam: numerical evaluation of the steady state differentiated w.r.t param_identif
%       nDdfstudt_Dparam: numerical evaluation of the degrees of freedom for student's t-distribution w.r.t param_identif
%
% Based upon Schmitt-Grohé and Uribe (2013) - Implementing Iskrev's Identifiability Test, http://www.columbia.edu/~mu2166/2nd_order.htm
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [ngra,nf,nhes,nSigma,netatilde,nsig,nSelectMat,nSS,ndfstudt,nd,nx,nv,ny,nu,...
    nDgra_Dparam,nDhes_Dparam,nDSigma_Dparam,nDetatilde_Dparam,nDsig_Dparam,nDsigetatilde_Dparam,nDSS_Dparam,nDdfstudt_Dparam]...
    = numeval(param_estim,DSGE_Model,approx,Deriv_type)

%% Evaluate numerically the steady state
eval(['[nSS] = ' DSGE_Model.shortname, '_SS(param_estim,DSGE_Model.param.names,DSGE_Model.spec);']);

%% Some precomputations
for i=1:length(nSS)
    eval(['vars' num2str(i) '= nSS(' num2str(i) ');']) % Set all variables equal to numerical steady-state
end
for i=1:length(param_estim)
    eval(['param' num2str(i) '= param_estim(' num2str(i) ');']) % Set all parameters equal to numerical local point
end

%% Evaluate NUMERICALLY the symbolic model objects in filename_num_eval.m
try
    eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(approx), '_num_eval']);
catch
    errordlg('Auxiliary Files have been saved to disk for the first time. Please run again! Sorry for that bug!');
end
%% Evaluate NUMERICALLY the analytical (symbolic) derivatives w.r.t. param_identif stored in filename_anal_deriv
if strcmp(Deriv_type,'Analytical')
    eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(approx), '_anal_deriv']);        
else
    % if numerical derivatives: don't calculate analytical derivatives
    nDgra_Dparam=[]; nDhes_Dparam=[]; nDSigma_Dparam=[]; nDetatilde_Dparam=[]; 
    nDsig_Dparam=[]; nDsigetatilde_Dparam=[]; nDSS_Dparam=[];  nDdfstudt_Dparam=[];   
end


end %main function end
