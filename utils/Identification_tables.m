function Identification_tables(DSGE_Model,Ident_Test,Settings)
% Tables with summary of rank test. Displays
%       - which Test, Model, Method, Order of Approximation, and Options for Test were used,
%       - which parameters were tested,
%       - which local point was used,
%       - checks order condition,
%       - Asks if the user wants to display problematic parameters using a singular value decomposition of the null space for the different
%       tolerance levels and objective matrices
%
% Inputs:
%       DSGE_Model: structure with information about DSGE model
%       Ident_Test: structure containing all relevant test results, i.e. ranks and problematic parameters
%       Settings: structure containing choosen Settings
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

% Some precomputations
table_results = Ident_Test.Results.tables.results;
steps = 9; tol0=1e-3;
table_bottom = Ident_Test.Results.tables.bottom;
pass_order = Ident_Test.Results.tables.pass_order;

%% Display results for rank tests
fprintf('\n\n');
fprintf('********** %s test **********\n',Ident_Test.fullname);
fprintf('Model: %s\n',DSGE_Model.fullname);
switch Settings.Derivative.type
    case 'Analytical'
        fprintf('Method: Analytical Derivatives\n');
    case 'Numerical'
        fprintf('Method: Numerical Derivatives, Differentiation Step: %s\n',Settings.Derivative.options{1});
end
fprintf('Order of Approximation: %d \n',Settings.approx);
fprintf('************************************************\n');
fprintf('Testing identification of:\n');
disp(DSGE_Model.param.identif_names);
fprintf('Local point:\n');
S = cell(2,size(DSGE_Model.param.estim,1));
for i=1:size(DSGE_Model.param.estim,1)
    S{1,i} = char(DSGE_Model.param.names(i));
    S{2,i} = DSGE_Model.param.estim(i);    
end
disp(S);
fprintf('**************** RESULT OF TEST ****************\n');

switch Ident_Test.shortname
    case 'Iskrev'
        fprintf(['Order condition: ' pass_order.Iskrev '\n']);
        fmt0=['%4e &' repmat('%5d&',1,size(table_results.Iskrev,2)-2) ' %5d' ' \\\\ \n'];
        fmt1=['%s&', repmat('%5d&',1, size(table_bottom.Iskrev,2)-1) ' %5d' ' \\\\ \n'] ;
        fprintf('  tol             J     M     Mbar  pass rank (Full Jacobian) \n' );
        for i=1:steps;
            fprintf(fmt0,table_results.Iskrev{i,:});
        end;
        fprintf(fmt1,' Robust      ', table_results.Iskrev{steps+1,2:end});
        fprintf(fmt1,' Require     ', table_bottom.Iskrev(1,:));
    case 'QuTkachenko'        
        fprintf(['Order condition: ' pass_order.QT '\n']);
        fmt0=['%4e &' repmat('%5d&',1,size(table_results.QT,2)-2) ' %5d' ' \\\\ \n'];
        fmt1=['%s&', repmat('%5d&',1, size(table_bottom.QT,2)-1) ' %5d' ' \\\\ \n'] ;
        fprintf(' tol              G   Gbar  pass rank \n' );
        for i=1:steps;
            fprintf(fmt0,table_results.QT{i,:});
        end;
        fprintf(fmt1,' Robust      ', table_results.QT{steps+1,2:end});
        fprintf(fmt1,' Require     ', table_bottom.QT(1,:));        
    case 'rank'
        fprintf(['Order condition Iskrev: ' pass_order.Iskrev '\n']);
        fprintf(['Order condition Qu & Tkachenko: ' pass_order.QT '\n']);
        fmt0=['%4e &' repmat('%5d&',1,size(table_results.Iskrev,2)-1+size(table_results.QT,2)-3) ' \\\\ \n'];
        fmt1=['%s&', repmat('%5d&',1, size(table_bottom.Iskrev,2)+size(table_bottom.QT,2)-2) ' \\\\ \n'] ;        
        fprintf(' tol              J     M     Mbar   G    Gbar \n' );
        for i=1:steps;
            fprintf(fmt0,[table_results.Iskrev{i,1:end-1} table_results.QT{i,2:end-1}]);
        end;    
        fprintf(fmt1,' Robust      ', [table_results.Iskrev{steps+1,2:end-1} table_results.QT{steps+1,2:end-1}]);
        fprintf(fmt1,' Require     ', [table_bottom.Iskrev(1,1:end-1) table_bottom.QT(1,1:end-1)]);
end

%% Display results for problematic parameters
print_probl_par = questdlg('Do you want to see the problematic parameters for each tolerance level',...
    'Problematic Parameters','Yes','No','Yes');
if strcmp(print_probl_par,'Yes')
    again = 'Yes';
    while strcmp(again,'Yes')
        % Which criteria to find bad?
        switch Ident_Test.shortname
            case 'Iskrev'
                choice = questdlg('For which case do you want to see the problematic parameters?', ...
                'Problematic Parameters', ...
                'Reduced-Form (J)','Cumulants Only (M)','Mean and Cumulants (Mbar)', 'Mean and Cumulants (Mbar)');
                switch choice
                    case 'Reduced-Form (J)'; tabl = Ident_Test.Results.ProblPars.Iskrev_J;
                    case 'Cumulants Only (M)'; tabl = Ident_Test.Results.ProblPars.Iskrev_M;
                    case 'Mean and Cumulants (Mbar)'; tabl = Ident_Test.Results.ProblPars.Iskrev_Mbar;
                end
            case 'QuTkachenko'
                choice = questdlg('For which case do you want to see the problematic parameters?', ...
                    'Problematic Parameters', ...
                        'Polyspectra only (G)','Mean and Polyspectra (Gbar)', 'Mean and Polyspectra (Gbar)');
                switch choice
                    case 'Polyspectra only (G)'; tabl = Ident_Test.Results.ProblPars.QT_G;
                    case 'Mean and Polyspectra (Gbar)'; tabl = Ident_Test.Results.ProblPars.QT_Gbar;
                end
            case 'rank'
                choice = bttnChoiseDialog({'Reduced-Form (J)' 'Cumulants Only (M)' 'Mean and Cumulants (Mbar)' 'Polyspectra only (G)' 'Mean and Polyspectra (Gbar)'}, 'Problematic Parameters', 'Mean and Polyspectra (Gbar)', 'For which case do you want to see the problematic parameters?');
                switch choice
                    case 1; tabl = Ident_Test.Results.ProblPars.Iskrev_J; choice='J';
                    case 2; tabl = Ident_Test.Results.ProblPars.Iskrev_M; choice='M';
                    case 3; tabl = Ident_Test.Results.ProblPars.Iskrev_Mbar; choice='Mbar';
                    case 4; tabl = Ident_Test.Results.ProblPars.QT_G; choice='G';
                    case 5; tabl = Ident_Test.Results.ProblPars.QT_Gbar; choice='Gbar';
                end
        end
        %% Display problematic parameters
        fprintf('** Problematic Parameters for %s **\n',choice);
        if isempty(tabl)
                fprintf('For all tolerance levels all parameters are distinguishable using %s\n',choice)                 
        else
            for i=1:(steps+1)
                % Set tolerance level
                if i <= steps
                    tol= tol0/10^(2*(i-1));
                else
                    tol = 'robust';
                end
                fprintf('***For tolerance level %s these parameter combinations are observationally equivalent using %s:\n\n',num2str(tol),choice)
                if i <= size(tabl,2)    
                    j=1; foundsome=0;
                    while j<=size(tabl{i},2)
                      if isempty(tabl{i}{j}) == 0
                          foundsome=1;
                          fprintf('%g-element subsets in each row:\n',j)
                          % Display subsets in each row. 
                          disp(tabl{i}{j}(:,:))                       
                      end              
                     j= j+1;
                    end
                end
                if foundsome ~= 1
                      fprintf('All parameters are distinguishable using %s \n',choice);
                end
                fprintf('\n');
            end
        end
        %Do you want other criteria?
        again = questdlg('Do you want to see another case','Problematic Parameters','Yes','No','No');        
    end
end

end %function end