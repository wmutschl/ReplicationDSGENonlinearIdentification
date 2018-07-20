% Get general settings in a GUI, i.e. which model, test and options, order of
% approximation, numerical or analytical derivatives
%
% Inputs: none
% Calls: bttnChoiseDialog
% Outputs:
%       DSGE_Model: Structure containing full name, short name and specification of model
%       Ident_Test: Structure containing options for identification tests
%       Settings: Structure containing Settings for analytical or numerical derivatives, approximation order
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [DSGE_Model,Ident_Test,Settings] = GUIGetSettings
% Some Options for GUI
options.Resize='on';                
options.Interpreter='tex';

% Choose model in a GUI (question dialog)
ChoiceModel = questdlg('Choose a model','Choose Model','An & Schorfheide (2007)','Kim (2003)','Kim (2003)');
switch ChoiceModel
    case 'An & Schorfheide (2007)'; DSGE_Model.shortname = 'AnSchorfheide'; DSGE_Model.fullname = 'An & Schorfheide (2007)';    
    case 'Kim (2003)'; DSGE_Model.shortname = 'Kim'; DSGE_Model.fullname = 'Kim (2003)';
end
addpath(['./models/', DSGE_Model.shortname],'-begin');

% Which specification of model
DSGE_Model.spec = str2double(inputdlg(sprintf(['Please choose the specification of the ' DSGE_Model.fullname ' Model: \n' eval([DSGE_Model.shortname '_spec'])]),'Specification',1,{'0'}));

% Choose order of approximation
ChoiceApprox = questdlg('Choose order of approximation','Approximation','1st order (not-efficient but okay)','2nd order','2nd order');
switch ChoiceApprox
    case '1st order (not-efficient but okay)'; Settings.approx=1;
    case '2nd order'; Settings.approx=2;
end

% Choose Identification procedure
Ident_Test.procedure.name = questdlg('Do you want to test a local point or prior domain?','Identification procedure','Local point','Prior domain','Local point');
if strcmp(Ident_Test.procedure.name,'Prior domain')
    Ident_Test.procedure.options = inputdlg({'Number of draws from prior','Rank tolerance (e.g. 1e-7 or robust)'},'Choose options for prior domain analysis',[1 60],{'100','robust'},options);        
end

% Choose Identification test
ChoiceTest = questdlg('Choose Identification Test','Identification Test','Both','Iskrev (2010)','Qu & Tkachenko (2012)','Both');

% Choose cumulants and polyspectra using checkboxes
cum.f = figure('units','normalized','position',[.3,0.5,.2,.1],...
             'toolbar','none','menu','none','Name','Which moments and spectra?','NumberTitle','off');
cum.second = uicontrol('Style','checkbox','units','normalized',...
                'String','Second moments/cumulants and power spectrum',...
                'Value',1,'Position',[0 0.8 1 0.2]);
cum.third = uicontrol('Style','checkbox','units','normalized',...
                'String','Third moments/cumulants and bispectrum',...
                'Value',1,'Position',[0 0.6 1 0.2]);
cum.fourth = uicontrol('Style','checkbox','units','normalized',...
                'String','Fourth moments/cumulants and trispectrum',...
                'Value',1,'Position',[0 0.4 1 0.2]);
cum.p = uicontrol('style','pushbutton','units','normalized',...% Create OK pushbutton
                'position',[0.4,0,.2,.2],'string','OK',...
                'callback',@p_call);
uiwait;

% Speed option: Compute cumulants of innovations and their derivatives only once
Settings.speed = questdlg('Select ''No Speed'' if you are running the model for (i) the first time, (ii) have changed the model variables or equations or (iii) plan to change the set of parameters to be identifiable compared to the last run. You need to run ''No Speed'' only once, after that you may choose ''Speed''!','Speed Option','Speed','No Speed','Speed');

% Choose method for problematic parameters
switch ChoiceTest
    case 'Both'
        Ident_Test.fullname = 'Criteria based on ranks';
        Ident_Test.shortname = 'rank';
        Ident_Test.options = inputdlg({'Method for problematic parameters: 1 for Brute Force (better) or 2 for Nullspace (faster)', 'Choose number of Lags in covariogram for Iskrev', 'Choose number of subintervalls for Qu & Tkachenko)'},'Options for rank tests',[1 90],{'1', '30' '100'},options);
    case 'Iskrev (2010)'
        Ident_Test.fullname = 'Iskrev (2010)';
        Ident_Test.shortname = 'Iskrev';
        Ident_Test.options = inputdlg({'Method for problematic parameters: 1 for Brute Force (better) or 2 for Nullspace (faster)','Choose number of Lags in covariogram'},'Options for Iskrev (2010)',[1 90],{'1' '30'},options);        
    case 'Qu & Tkachenko (2012)'
        Ident_Test.fullname = 'Qu & Tkachenko (2012)';
        Ident_Test.shortname = 'QuTkachenko';
        Ident_Test.options = inputdlg({'Method for problematic parameters: 1 for Brute Force (better) or 2 for Nullspace (faster)','Choose number of subintervalls'},'Options for Qu & Tkachenko (2012)',[1 90],{'1' '100'},options);
end

% Choose type of derivatives
Settings.Derivative.type = questdlg({'Do you want analytical or numerical derivatives?'},'Derivatives','Analytical','Numerical','Analytical');
if strcmp(Settings.Derivative.type,'Numerical')   	
    Settings.Derivative.options = inputdlg({'Choose differentiation step'},'Options for numerical derivatives',[1 50],{num2str(1e-7)},options);
end

% Set folders for auxiliary model files, dependent on the choosen model
addpath(['./models/', DSGE_Model.shortname],'-begin');


function p_call(varargin)
    % Get values for checkboxes    
    Ident_Test.cumulants.second = get(cum.second,'Value'); 
    Ident_Test.cumulants.third = get(cum.third,'Value'); 
    Ident_Test.cumulants.fourth = get(cum.fourth,'Value'); 
    close(gcf)
end
end %main function end






