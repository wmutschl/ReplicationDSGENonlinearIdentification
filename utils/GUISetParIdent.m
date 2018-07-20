% Select for which parameters identifiabilty is to be checked, make them independent of
% names, set local point. All in a Graphical User Interface(GUI).
%
% Inputs: 
%       params: structure containing all information about parameters of choosen model in original form
% Outputs: 
%       param: structure containing all information about parameters of choosen model in independent-of-names form
% Calls:
%       uibutton, plot_priors
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function param = GUISetParIdent(params)
% Initialize variables
nparam =size(params,2); 
param.names =[]; param.names_tex={};
param.default.point={}; param.prior.type = {};param.prior.par1 = {};
param.prior.par2 = {}; param.bound.type = {}; param.bound.low = {}; 
param.bound.up =   {}; param.bound.c = {}; 

% Get values from structure in original form
for i = 1:nparam   
    param.names =        [param.names params{i}{1}];
    param.names_tex=     [param.names_tex params{i}{2}];
    param.default.point= [param.default.point params{i}{3}];
    param.prior.type =   [param.prior.type params{i}{4}];
    param.prior.par1 =   [param.prior.par1 params{i}{5}];
    param.prior.par2 =   [param.prior.par2 params{i}{6}];
    param.bound.type =   [param.bound.type params{i}{7}];
    param.bound.low =    [param.bound.low params{i}{8}];
    param.bound.up =     [param.bound.up params{i}{9}];
    param.bound.c =      [param.bound.c params{i}{10}];
end

% Create strings for original parameter names
param_names_string = cell(size(param.names));
for i=1:length(param_names_string)
    param_names_string{i} = char(param.names(i));
end    


%% Graphical User Interface
% Set local point, i.e. assign a value to all parameters of the DSGE model (not just those whose identifiability is being checked)
% Check the parameters whose identifiability is being checked
     % For the figure only 18 parameters are displayed in one column
     if (nparam/18-floor(nparam/18)) == 0
         ncols  = nparam/18;
     else
         ncols = floor(nparam/18)+1;
     end     
   
    % Create figure, you can change its size, this works for up to 36 parameters
    h.f = figure('units','normalized','position',[.3,0.2,.3,.75],...
             'toolbar','none','menu','none','Name','Set parameters (all values may be changed here!)','NumberTitle','off');
    
    x = nparam; % initialize auxiliary index for the loop for each column
    for ncol=1:ncols
        % Create titles
        uibutton('style','text','units','normalized',...
                'position',[0.02+(ncol-1)*.5,0.95,.12,.03],'String','Parameter');
        uibutton('style','text','units','normalized',...
                'position',[0.15+(ncol-1)*.5,0.95,.12,.03],'String','Local Point');
        uibutton('style','text','units','normalized',...
                'position',[0.28+(ncol-1)*.5,0.95,.12,.03],'String','Prior Type');
        uibutton('style','text','units','normalized',...
                'position',[0.41+(ncol-1)*.5,0.95,.12,.03],'String','Prior Par 1');        
        uibutton('style','text','units','normalized',...
                'position',[0.54+(ncol-1)*.5,0.95,.12,.03],'String','Prior Par 2');        
        uibutton('style','text','units','normalized',...
                'position',[0.67+(ncol-1)*.5,0.95,.12,.03],'String','Lower Bound');
        uibutton('style','text','units','normalized',...
                'position',[0.80+(ncol-1)*.5,0.95,.12,.03],'String','Upper Bound');        
        uibutton('style','text','units','normalized',...
                'position',[0.93+(ncol-1)*.5,0.95,.12,.03],'String','Check?');        

        % Settings for loop if more than one column
        if x > 18
            iend = ncol*18;
            x = x-18;
        else
            iend = nparam;
        end
        j=1; % auxiliary index for the rows in each column
        for i=(1+(ncol-1)*18):iend
            h.parnam(i) = uibutton('style','text','units','normalized',...
                    'position',[0.02+(ncol-1)*.5,0.95-j*0.05,.12,.03],'String',param.names_tex(i),'Interpreter','tex');
            h.localpoint(i) = uicontrol('style','edit','units','normalized',...
                    'position',[0.15+(ncol-1)*.5,0.95-j*0.05,.12,.03],'String',param.default.point(i));            
            h.priortype(i) = uicontrol('style','edit','units','normalized',...
                    'position',[0.28+(ncol-1)*.5,0.95-j*0.05,.12,.03],'String',param.prior.type(i));
            h.priorpar1(i) = uicontrol('style','edit','units','normalized',...
                    'position',[0.41+(ncol-1)*.5,0.95-j*0.05,.12,.03],'String',param.prior.par1(i));
            h.priorpar2(i) = uicontrol('style','edit','units','normalized',...
                    'position',[0.54+(ncol-1)*.5,0.95-j*0.05,.12,.03],'String',param.prior.par2(i));
            h.priorboundlow(i) = uicontrol('style','edit','units','normalized',...
                    'position',[0.67+(ncol-1)*.5,0.95-j*0.05,.12,.03],'String',param.bound.low(i));
            h.priorboundup(i) = uicontrol('style','edit','units','normalized',...
                    'position',[0.80+(ncol-1)*.5,0.95-j*0.05,.12,.03],'String',param.bound.up(i));
            h.check(i) = uicontrol('style','checkbox','units','normalized',...
                    'position',[0.93+(ncol-1)*.5,0.95-j*0.05,.12,.03],'Value',1);
            j=j+1;
        end
    end
    % Create OK pushbutton
    h.p = uicontrol('style','pushbutton','units','normalized',...
                'position',[.5,0,.2,.05],'string','OK',...
                'callback',@p_call);

    uiwait;

%% Make identifiable parameters independent of names and save local point in a vector
param.identif = sym(zeros(length(param.identif_names),1));
param.estim = zeros(nparam,1);
param.fix = ones(nparam,1);
param.numbered=sym(zeros(1,length(param.identif)));
for i=1:nparam
    param.numbered(i) = sym(['param' num2str(i)]);
    param.estim(i) = str2double(localpoint{i});
    param.prior.type(i) = priortype{i};
    param.prior.par1{i} = str2double(priorpar1{i});
    param.prior.par2{i} = str2double(priorpar2{i});
    param.bound.low{i} = str2double(priorboundlow{i}); 
    param.bound.up{i} = str2double(priorboundup{i}); 
    for j=1:length(param.identif_names)
        if sym(param.identif_names(j)) == param.names(i)
            param.identif(j) = ['param' num2str(i)];
            param.fix(i) = 0;
        end        
    end
end

%% Plot priors
if strcmp(questdlg('Do you want to plot the priors?','Plotting Priors','Yes','No','No'),'Yes')
    plot_priors(param);
end

%% Auxiliary functions for GUI: Pushbutton callback
function p_call(varargin)    
    vals = get(h.check,'Value'); 
    param.identif_names=param_names_string(find([vals{:}])); 
    localpoint =get(h.localpoint,'String');
    priortype = get(h.priortype,'String');
    priorpar1 = get(h.priorpar1,'String');
    priorpar2 = get(h.priorpar2,'String');
    priorboundlow = get(h.priorboundlow,'String'); 
    priorboundup =  get(h.priorboundup,'String'); 
    close(gcf)
end

end
