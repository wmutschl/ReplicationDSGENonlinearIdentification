%% Matlab Code for "Identification of DSGE Models - the Effect of Higher-Order Approximation and Pruning"
%
% Compares identification tests for two example models:
% 'An & Schorfheide (2007)' 'Kim (2003)'
%
% User can choose in a graphical-user-interface between the models, the tests, which parameters to
% check at which local point or specifying a prior domain, analytical or numerical derivatives, 
% the order of approximation (up to second), which cumulants/moments (up to
% fourth) and corresponding polyspectra to consider.
% Note: Order=1 is not efficient, since it uses the pruned state-space representation as well
%
% Since all procedures are model independent, other models can be easily
% included and tested as long as they are put into the same framework.
% 
% Inputs: None
%
% Output: Table or Graphs with identification results and analysis of problematic parameters
%
% Just run identfication_run.m.
% Note: Matlab's symbolic toolbox is needed. 
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

%%
% Housekeeping and some options
clear all
addpath('./utils','./models','-begin');   
% Get general Settings in a GUI i.e. which model, test and options, order of approximation, numerical or analytical derivatives
[DSGE_Model,Ident_Test,Settings] = GUIGetSettings;
% GUI: Select which parameters are to be identifiable, make them independent of names, set local point
DSGE_Model.param = GUISetParIdent(eval([DSGE_Model.shortname '_param(DSGE_Model.spec);']));

%% Preprocessing model: put into desired DSGE framework independent of names (parameters and variables are numbered)
% Get symbolic expressions for model variables and objects all evaluated at the steady state, print all expressions to file 
[DSGE_Model.symbolic,DSGE_Model.numbers] = ModelPostCalculations(eval([DSGE_Model.shortname '(DSGE_Model.spec)']),DSGE_Model.param,Settings);
% Check existence of auxiliary files to prevent bug if you run the model for the first time
CheckFileExistence(DSGE_Model,Settings.approx);
% Print and save all expressions to file `filename_num_eval.m' for further evaluation
anal_deriv_f_print2f(DSGE_Model,Settings.approx);

% Get symbolic expression for steady-state
fprintf('GET SYMBOLIC EXPRESSION FOR STEADY STATE...');
DSGE_Model.symbolic.SS = eval([DSGE_Model.shortname '_SS(DSGE_Model.param.numbered,DSGE_Model.param.names,DSGE_Model.spec);']);
fprintf('FINISHED!\n');
% Prints to a file the symbolic derivatives of model objects w.r.t. param_identif for further evaluation
if strcmp(Settings.Derivative.type,'Analytical')
    anal_deriv_param(DSGE_Model,Settings.approx);
end

%% Run choosen identification test and print results
fprintf('Start checking %s criteria \n',Ident_Test.fullname);
Ident_Test.Results = rank_tests(DSGE_Model,Settings,Ident_Test);

%% Display results of different identification tests and problematic parameters
switch Ident_Test.procedure.name
    case 'Prior domain'
        Identification_barplots(Ident_Test);
    case 'Local point'
        Identification_tables(DSGE_Model,Ident_Test,Settings);
end

%% Clean up
rmpath('./utils','./models',['./models/', DSGE_Model.shortname]);
 