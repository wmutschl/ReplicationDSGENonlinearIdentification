% Displays barplots for problematic parameters of identification tests of the prior domain
%
% Inputs:       Ident_Test: Structure containing results of tests for the prior domain
%
% Calls:        count_unique
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function Identification_barplots(Ident_Test)
num_draws= str2double(Ident_Test.procedure.options{1});

switch Ident_Test.shortname
    case 'Iskrev'
        PROBL_PARS = {Ident_Test.Results.ProblPars.Iskrev_J Ident_Test.Results.ProblPars.Iskrev_M Ident_Test.Results.ProblPars.Iskrev_Mbar};
        TITLES = {'J' 'M' 'Mbar'};
    case 'QuTkachenko'
        PROBL_PARS = {Ident_Test.Results.ProblPars.QT_G Ident_Test.Results.ProblPars.QT_Gbar};
        TITLES = {'G' 'Gbar'};
    case 'rank'
        PROBL_PARS = {Ident_Test.Results.ProblPars.Iskrev_J Ident_Test.Results.ProblPars.Iskrev_M Ident_Test.Results.ProblPars.Iskrev_Mbar...
              Ident_Test.Results.ProblPars.QT_G Ident_Test.Results.ProblPars.QT_Gbar};
        TITLES = {'J' 'M' 'Mbar' 'G' 'Gbar'};
end

for k = 1:size(PROBL_PARS,2)
    probl_pars=[];
    if ~isempty(PROBL_PARS{k}(~cellfun('isempty',PROBL_PARS{k})))
        fprintf('Sets responsible for rank failure in percentage of total draws using %s:',TITLES{k});
        for num_draw = 1:num_draws
            for num_set=1:size(PROBL_PARS{k}{num_draw},2) %For j-element subsets
                S={};
                for i=1:size(PROBL_PARS{k}{num_draw}{num_set},1)
                    S{i} = strjoin(PROBL_PARS{k}{num_draw}{num_set}(i,:),',');
                end
                probl_pars = [probl_pars S];
            end
        end
        [uniques,numUnique] = count_unique(probl_pars);        
        table(uniques,numUnique/num_draws*100,'VariableNames',{'Sets' 'Percentage_of_draws_that_are_not_identified'})
        X = [numUnique'/num_draws; numUnique'*0]*100;   % trick: add a second category of all zeros
        figure
        bar(X);         % now we have a plot with an ugly second category
        xlim([0.5 1.5]); % trick that hides second category of X
        ylim([0 100]);
        legend(uniques,'Location','BestOutside'); % add the legend
        set(gca, 'XTickLabel', ''); % hide the '1' on the X axis (or make it whatever you want)    
        title(['Sets responsible for rank failure in percentage of total draws using ' TITLES{k}]);
    else
        fprintf('All parameters are identified, using the %s criteria!\n',TITLES{k});
    end
end
