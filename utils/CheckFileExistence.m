% Checks if auxiliary script files are existant. If not it creates them empty, otherwise it gives a strange bug for the first time a model is computed.
%
% Inputs:
%       DSGE_Model: structure with information about DSGE model
%       approx: Order of approximation
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)
function CheckFileExistence(DSGE_Model,approx)
spec = DSGE_Model.spec;
model = DSGE_Model.shortname;
filename = ['./models/', DSGE_Model.shortname,'/',DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx',];
if ~exist([filename,num2str(approx),'_num_eval.m'], 'file')
    fid1=fopen([filename,num2str(approx),'_num_eval.m'],'w');
    fclose(fid1);
end

if ~exist(['./models/',model,'/', model,'_spec',num2str(spec),'_approx',num2str(approx),'_anal_deriv.m'],'file')
    fid2=fopen(['./models/',model,'/', model,'_spec',num2str(spec),'_approx',num2str(approx),'_anal_deriv.m'],'w');
    fclose(fid2);
end

filename = ['./models/', DSGE_Model.shortname,'/',DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx',];

if ~exist([filename,num2str(approx),'_prodmom2_num_eval.m'], 'file')
    fid3=fopen([filename,num2str(approx),'_prodmom2_num_eval.m'],'w');
    fclose(fid3);
end
if ~exist([filename,num2str(approx),'_prodmom3_num_eval.m'], 'file')
    fid4=fopen([filename,num2str(approx),'_prodmom3_num_eval.m'],'w');
    fclose(fid4);
end
if ~exist([filename,num2str(approx),'_prodmom4_num_eval.m'], 'file')
    fid5=fopen([filename,num2str(approx),'_prodmom4_num_eval.m'],'w');
    fclose(fid5);
end
if ~exist([filename,num2str(approx),'_prodmom2_deriv.m'], 'file')
    fid6=fopen([filename,num2str(approx),'_prodmom2_deriv.m'],'w');
    fclose(fid6);
end
if ~exist([filename,num2str(approx),'_prodmom3_deriv.m'], 'file')
    fid7=fopen([filename,num2str(approx),'_prodmom3_deriv.m'],'w');
    fclose(fid7);
end
if ~exist([filename,num2str(approx),'_prodmom4_deriv.m'], 'file')
    fid8=fopen([filename,num2str(approx),'_prodmom4_deriv.m'],'w');
    fclose(fid8);
end

fclose('all');   
end

