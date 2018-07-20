% Function that computes identification tests based on ranks (Iskrev, Qu&Tkachenko) for
% different tolerance levels and analytical or numerical derivatives, and first and second order approximation.
% Note: Procedures use pruned state-space also for first-order approximation (not efficient)!
% Finds also problematic parameters if ranks are deficient.
% Procedures are either for local point or prior domain.
%
% Input:
%       DSGE_Model:structure containing all information about DSGE model independent of names
%       Settings: structure containing Settings for derivatives and order of approximation
%       Ident_Test: structure containing settings for identification tests
% Calls:
%       drawprior, AnalyticalDerivatives, NumericalDerivatives, Nulltol,
%       Find_Problematic_Params
% Output:
%       Test_Results: structure containing Jacobians and rank results of tests
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function Test_Results = rank_tests(DSGE_Model,Settings,Ident_Test)
proc = Ident_Test.options{1}; % Settings for procedure to find problematic parameters
parnames = transpose(DSGE_Model.param.identif_names); % names of parameters
nparam = length(parnames); % number of parameters
    
if strcmp(Ident_Test.procedure.name, 'Prior domain')
    num_draws = str2double(Ident_Test.procedure.options{1}); % number of draws from prior domain
    steps=0; % no steps in rank calculation, i.e. only one tolerance level
    if strcmp(Ident_Test.procedure.options{2},'robust')
        tol_spec = Ident_Test.procedure.options{2};
    else
        tol_spec = str2double(Ident_Test.procedure.options{2});
    end 
else
    num_draws = 1; tol0=1e-3; steps=9; % local point: only one draw, different tolerance level starting with tol0 and steps
end

%% Initialize storage
%Iskrev
if xor(strcmp(Ident_Test.shortname,'Iskrev'),strcmp(Ident_Test.shortname,'rank'))
    ProblPars.Iskrev_J=cell(num_draws,steps+1);
    ProblPars.Iskrev_M=cell(num_draws,steps+1);
    ProblPars.Iskrev_Mbar=cell(num_draws,steps+1);
    T = str2double(Ident_Test.options{2}); %Horizon ob lags in Iskrev
end
%Qu and Tkachenko
if xor(strcmp(Ident_Test.shortname,'QuTkachenko'),strcmp(Ident_Test.shortname,'rank'))
    ProblPars.QT_G=cell(num_draws,steps+1);
    ProblPars.QT_Gbar=cell(num_draws,steps+1);
end

for num_draw = 1:num_draws
    %% Set local point
    if strcmp(Ident_Test.procedure.name, 'Prior domain')
        fprintf('Checking for prior domain. Starting %d out of %d\n',num_draw,num_draws);
        retcode =0; 
        while retcode == 0
            [param_estim,retcode] = drawprior(DSGE_Model,Settings);            
            if retcode==0 % Use only draws from determinacy region
                fprintf('...discard prior draw...get new one...\n');
            end
        end
    else
        fprintf('Checking local point...\n');
        param_estim = DSGE_Model.param.estim;      % Set local point as specified in the GUI  
    end
   
    %% Compute derivatives either analytically or numerically
    % Analytical Derivatives
    if strcmp(Settings.Derivative.type,'Analytical')        
        fprintf('    Calculating Solution and Analytical Derivatives...\n');
        [Solut,Deriv] = AnalyticalDerivatives(param_estim,DSGE_Model,Settings,Ident_Test);
        fprintf('Calculating Solution and Analytical Derivatives finished!\n');
    % Numerical Derivatives
    elseif strcmp(Settings.Derivative.type,'Numerical')
        fprintf('    Calculating Solution and Numerical Derivatives...\n');
        [Solut,Deriv] = NumericalDerivatives(param_estim,DSGE_Model,Settings,Ident_Test);
        fprintf('Calculating Solution and Numerical Derivatives finished!\n');
    end
    
    
    %% Set up Iskrev 
    % Iskrev Jacobians
    if xor(strcmp(Ident_Test.shortname,'Iskrev'),strcmp(Ident_Test.shortname,'rank'))
        fprintf('    Calculating Jacobians for Iskrev...');
        Jacobians.Iskrev_J = full([...
             Solut.SelectMat*Deriv.DSS_Dparam(DSGE_Model.numbers.nv+1:end,:);
             Deriv.Dprun_c_Dparam;
             Deriv.Dprun_d_Dparam;
             Deriv.Dprun_A_Dparam;
             Deriv.Dprun_B_Dparam;
             Deriv.Dprun_C_Dparam;
             Deriv.Dprun_D_Dparam;
             Deriv.DGAMMA2min_Dparam;
             Deriv.DGAMMA3min_Dparam;
             Deriv.DGAMMA4min_Dparam]);
        Jacobians.Iskrev_M = full([Deriv.Iskrev_M]);
        Jacobians.Iskrev_Mbar = full([Deriv.DEd_Dparam;Deriv.Iskrev_M]);

        %Rescale rows by its largest element in absolute value, rank stays the same        
        indx_J = find(sum(Jacobians.Iskrev_J,2)~=0);
        indx_M = find(sum(Jacobians.Iskrev_M,2)~=0);
        indx_Mbar = find(sum(Jacobians.Iskrev_Mbar,2)~=0);
        Jacobians.Iskrev_J(indx_J,:) = Jacobians.Iskrev_J(indx_J,:)./repmat(max(abs(Jacobians.Iskrev_J(indx_J,:)'))',1,size(Jacobians.Iskrev_J,2));
        Jacobians.Iskrev_M(indx_M,:) = Jacobians.Iskrev_M(indx_M,:)./repmat(max(abs(Jacobians.Iskrev_M(indx_M,:)'))',1,size(Jacobians.Iskrev_M,2));
        Jacobians.Iskrev_Mbar(indx_Mbar,:) = Jacobians.Iskrev_Mbar(indx_Mbar,:)./repmat(max(abs(Jacobians.Iskrev_Mbar(indx_Mbar,:)'))',1,size(Jacobians.Iskrev_Mbar,2));
    
        % Iskrev tables
        n_required_order_Iskrev= (T-1)*DSGE_Model.numbers.nd^2+DSGE_Model.numbers.nd*(DSGE_Model.numbers.nd+3)/2;
        if (nparam <= n_required_order_Iskrev) ==1
            tables.pass_order.Iskrev = 'pass';
        else
            tables.pass_order.Iskrev = 'fail'; 
        end
        tables.bottom.Iskrev = [nparam nparam nparam  1];
        fprintf('finished!\n');
    end
    
    %% Set up Qu and Tkachenko 
    % Qu and Tkachenko Jacobians
    if xor(strcmp(Ident_Test.shortname,'QuTkachenko'),strcmp(Ident_Test.shortname,'rank'))
        fprintf('    Calculating Jacobians for Qu and Tkachenko...');        
        Jacobians.QT_G = full(Deriv.QT_G);
        Jacobians.QT_Gbar = Deriv.QT_G + full(transpose(Deriv.DEd_Dparam)*Deriv.DEd_Dparam);
        %Rescale rows by its largest element in absolute value, order of magnitude stays the same
        indx_G = find(sum(Jacobians.QT_G,2)~=0);
        indx_Gbar = find(sum(Jacobians.QT_Gbar,2)~=0);    
        Jacobians.QT_G(indx_G,:) = Jacobians.QT_G(indx_G,:)./repmat(max(abs(Jacobians.QT_G(indx_G,:)'))',1,size(Jacobians.QT_G,2));
        Jacobians.QT_Gbar(indx_Gbar,:) = Jacobians.QT_Gbar(indx_Gbar,:)./repmat(max(abs(Jacobians.QT_Gbar(indx_Gbar,:)'))',1,size(Jacobians.QT_Gbar,2));
        % Qu and Tkachenko tables
        tables.bottom.QT=[nparam nparam 1];
        tables.pass_order.QT = 'No order condition';
        fprintf('finished!\n');
    end
    
    %% Calculating ranks and problematic parameters for different tolerance levels
    for i=1:(steps+1)
        % Set tolerance level
        if strcmp(Ident_Test.procedure.name, 'Local point')
            if i <= steps
                tol= tol0/10^(2*(i-1));
                fprintf('Calculating ranks and problematic parameters for tolerance level %e\n',tol)
            else
                tol = 'robust'; %Last tolerance level in loop is robust level for local point procedure
                fprintf('Calculating ranks and problematic parameters for tolerance level %s\n',tol)
            end
        elseif strcmp(Ident_Test.procedure.name, 'Prior domain')
            tol = tol_spec;
            fprintf('    Calculating ranks and problematic parameters for tolerance level %s:\n',tol)
        end    
    
        % Iskrev ranks and problematic parameters
        if xor(strcmp(Ident_Test.shortname,'Iskrev'),strcmp(Ident_Test.shortname,'rank'))
            fprintf('    Iskrev...');
            % Compute rank and check failures of J
            [n_J,Null_J,tol_J] = NullTol(Jacobians.Iskrev_J,tol);
            n_slack_J = nparam-n_J;            
            if n_slack_J> 0
                ProblPars.Iskrev_J{num_draw,i} = Find_Problematic_Params(Jacobians.Iskrev_J,Null_J,tol_J,parnames,nparam,proc,'Iskrev');
            end  

            % Compute rank and check failures of Moments only, i.e. M
            [n_M,Null_M,tol_M] = NullTol(Jacobians.Iskrev_M,tol);        
            n_slack_M = nparam-n_M;            
            if n_slack_M> 0
                ProblPars.Iskrev_M{num_draw,i} = Find_Problematic_Params(Jacobians.Iskrev_M,Null_M,tol_M,parnames,nparam,proc,'Iskrev');
            end

            % Compute rank and check failures of mean and moments, i.e.Mbar
            [n_Mbar,Null_Mbar,tol_Mbar] = NullTol(Jacobians.Iskrev_Mbar,tol);
            n_slack_Mbar = nparam-n_Mbar;
            if n_slack_Mbar> 0
                ProblPars.Iskrev_Mbar{num_draw,i} = Find_Problematic_Params(Jacobians.Iskrev_Mbar,Null_Mbar,tol_Mbar,parnames,nparam,proc,'Iskrev');
            end

            % Check rank condition for Mbar
            pass_rank_Iskrev = (n_Mbar == nparam);
            % Save ranks in table
            tables.results.Iskrev(i,:)={tol n_J  n_M  n_Mbar pass_rank_Iskrev};
        end
        
        % Qu and Tkachenko ranks and problematic parameters
        if xor(strcmp(Ident_Test.shortname,'QuTkachenko'),strcmp(Ident_Test.shortname,'rank'))
            fprintf('    Qu & Tkachenko...');
            % Compute rank of Gbar
            [n_Gbar,Null_Gbar,tol_Gbar] = NullTol(Jacobians.QT_Gbar,tol);
            n_slack_Gbar = nparam-n_Gbar;            
            if n_slack_Gbar> 0            
                ProblPars.QT_Gbar{num_draw,i} = Find_Problematic_Params(Jacobians.QT_Gbar,Null_Gbar,tol_Gbar,parnames,nparam,proc,'QuTkachenko');
            end

            % Compute rank of G
            [n_G,Null_G,tol_G] = NullTol(Jacobians.QT_G,tol);
            n_slack_G = nparam-n_G;            
            if n_slack_G> 0
                ProblPars.QT_G{num_draw,i} = Find_Problematic_Params(Jacobians.QT_G,Null_G,tol_G,parnames,nparam,proc,'QuTkachenko');        
            end

            % Check rank condition
            pass_rank_QT = (n_Gbar == nparam);
            % Save ranks in table
            tables.results.QT(i,:)={tol n_G  n_Gbar pass_rank_QT};
        end

        fprintf('finished!\n');
    end %steps end

end %numdraws end

%Store Results of Test
Test_Results.Jacobians = Jacobians;
Test_Results.tables =tables;
Test_Results.ProblPars = ProblPars;



end %main function end