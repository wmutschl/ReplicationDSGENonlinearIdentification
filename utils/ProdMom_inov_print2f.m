% This function writes all product moments and analytical derivatives of
% xi_t into files
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function ProdMom_inov_print2f(auxpar,DSGE_Model,approx,Deriv_type,Ident_Test,M2,M3,M4,DM2_Dauxpar,DM3_Dauxpar,DM4_Dauxpar,ic2,ic3,ic4)
filename = ['./models/', DSGE_Model.shortname,'/',DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx',];
funcname = [DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx'];
for i=1:length(auxpar)        
    arg{i} = char(auxpar(i));
end    

%% Second-order Product Moments
if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
fprintf('WRITE SYMBOLIC EXPRESSIONS OF SECOND-ORDER PRODUCT MOMENTS OF INNOVATIONS TO M-FILE...');
% Save into file
fid=fopen([filename,num2str(approx),'_prodmom2_num_eval.m'],'w');
fprintf(fid,'function [nM2,ic2] = %s%d_prodmom2_num_eval(arg) \n',funcname,approx);
for i=1:length(arg)
    fprintf(fid,'%s = arg(%d); \n',arg{i},i);            
end
% PRINT M2
    izero=(M2~=0);
    izeroi=find(izero);
    S=['nM2=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi)
     S2 = char(M2(izeroi(i)));     
     S = ['nM2(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nM2=reshape(nM2,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
    
% Print ic2
    S = mat2str(ic2);
    S=['ic2=', S,';\n'];
    S=regexprep(S, ',', '\r');    
    fprintf(fid,S);    

fclose(fid);
fprintf('FINISHED!\n');
end

%% Third-order Product Moments
if Ident_Test.cumulants.third
fprintf('WRITE SYMBOLIC EXPRESSIONS OF THIRD-ORDER PRODUCT MOMENTS OF INNOVATIONS TO M-FILE...');
% Save into file
fid=fopen([filename,num2str(approx),'_prodmom3_num_eval.m'],'w');
fprintf(fid,'function [nM3,ic3] = %s%d_prodmom3_num_eval(arg) \n',funcname,approx);
for i=1:length(arg)
    fprintf(fid,'%s = arg(%d); \n',arg{i},i);            
end
% PRINT M3
    izero=(M3~=0);
    izeroi=find(izero);
    S=['nM3=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi)
     S2 = char(M3(izeroi(i)));
     S = [' nM3(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nM3=reshape(nM3,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 

% Print ic3
    S = mat2str(ic3);
    S=['ic3=', S,';\n'];
    S=regexprep(S, ',', '\r');    
    fprintf(fid,S);    

fclose(fid);
fprintf('FINISHED!\n');
end

%% Fourth-order Product Moments
if Ident_Test.cumulants.fourth
fprintf('WRITE SYMBOLIC EXPRESSIONS OF FOURTH-ORDER PRODUCT MOMENTS OF INNOVATIONS TO M-FILE...');
% Save into file
fid=fopen([filename,num2str(approx),'_prodmom4_num_eval.m'],'w');
fprintf(fid,'function [nM4,ic4] = %s%d_prodmom4_num_eval(arg) \n',funcname,approx);
for i=1:length(arg)
    fprintf(fid,'%s = arg(%d); \n',arg{i},i);            
end
   
% PRINT M4
    izero=(M4~=0);
    izeroi=find(izero);
    S=['nM4=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi)
     S2 = char(M4(izeroi(i)));
     S = [' nM4(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nM4=reshape(nM4,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 

% Print ic4
    S = mat2str(ic4);
    S=['ic4=', S,';\n'];
    S=regexprep(S, ',', '\r');    
    fprintf(fid,S);    

fclose(fid);
fprintf('FINISHED!\n');
end

if strcmp(Deriv_type,'Analytical')
%% ANALYTICAL DERIVATIVES OF SECOND-ORDER PRODUCT MOMENTS
if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
fprintf('WRITE SYMBOLIC EXPRESSIONS OF DERIVATIVES OF SECOND-ORDER PRODUCT MOMENTS OF INNOVATIONS TO M-FILE...')
fid=fopen([filename,num2str(approx),'_prodmom2_deriv.m'],'w');
fprintf(fid,'function nDM2_Dauxpar = %s%d_prodmom2_deriv(arg) \n',funcname,approx);
for i=1:length(arg)
    fprintf(fid,'%s = arg(%d); \n',arg{i},i);            
end
    % Print DM2_Dauxpar
    izero=(DM2_Dauxpar~=0);
    izeroi=find(izero);
    S=['nDM2_Dauxpar=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi)
     S2 = char(DM2_Dauxpar(izeroi(i)));
     S = [' nDM2_Dauxpar(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDM2_Dauxpar=reshape(nDM2_Dauxpar,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
fclose(fid);
fprintf('FINISHED!\n');
end

%% ANALYTICAL DERIVATIVES OF THIRD-ORDER PRODUCT MOMENTS
if Ident_Test.cumulants.third
fprintf('WRITE SYMBOLIC EXPRESSIONS OF DERIVATIVES OF THIRD-ORDER PRODUCT MOMENTS OF INNOVATIONS TO M-FILE...')
fid=fopen([filename,num2str(approx),'_prodmom3_deriv.m'],'w');
fprintf(fid,'function nDM3_Dauxpar = %s%d_prodmom3_deriv(arg) \n',funcname,approx);
for i=1:length(arg)
    fprintf(fid,'%s = arg(%d); \n',arg{i},i);            
end
    izero=(DM3_Dauxpar~=0);
    izeroi=find(izero);
    S=['nDM3_Dauxpar=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi)
     S2 = char(DM3_Dauxpar(izeroi(i)));
     S = [' nDM3_Dauxpar(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDM3_Dauxpar=reshape(nDM3_Dauxpar,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
fclose(fid);
fprintf('FINISHED!\n');
end

%% ANALYTICAL DERIVATIVES OF FOURTH-ORDER PRODUCT MOMENTS
if Ident_Test.cumulants.fourth
fprintf('WRITE SYMBOLIC EXPRESSIONS OF DERIVATIVES OF FOURTH-ORDER PRODUCT MOMENTS OF INNOVATIONS TO M-FILE...')
fid=fopen([filename,num2str(approx),'_prodmom4_deriv.m'],'w');
fprintf(fid,'function nDM4_Dauxpar = %s%d_prodmom4_deriv(arg) \n',funcname,approx);
for i=1:length(arg)
    fprintf(fid,'%s = arg(%d); \n',arg{i},i);            
end
    izero=(DM4_Dauxpar~=0);
    izeroi=find(izero);
    S=['nDM4_Dauxpar=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi)
     S2 = char(DM4_Dauxpar(izeroi(i)));
     S = [' nDM4_Dauxpar(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDM4_Dauxpar=reshape(nDM4_Dauxpar,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 

fclose(fid);
fprintf('FINISHED!\n');
end
end %analytical derivatives end
