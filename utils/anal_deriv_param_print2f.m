% Writes symbolic derivatives of model objects w.r.t. param_identif to an m-file for numeric evaluation
%
% Inputs: 
%       model: short name of model used for filename
%       spec: specification of model
%       approx: order of approximation
%       DSS_Dparam: analytical derivative of steady-state w.r.t. param_identif
%       Df_Donlyparam: analytical derivative of w.r.t. param_identif only
%       Df_Dvars: analytical derivative of f w.r.t. steady state vector
%       Dgra_Donlyparam: analytical derivative of gradient of f w.r.t. param_identif only
%       Dgra_Dvars: analytical derivative of gradient of f w.r.t. steady state vector
%       Dhes_Donlyparam: analytical derivative of Magnus-Neudecker-Hessian of f w.r.t. param_identif only
%       Dhes_Dvars: analytical derivative of Magnus-Neudecker-Hessian of f w.r.t. steady state vector
%       DSigma_Dparam: analytical derivative of var/cov matrix of shocks and measurement errors in transition equation including perturbation parameter w.r.t. param_identif
%       Dsig_Dparam: analytical derivative of perturbation paramter w.r.t. param_identif
%       Dsigetatilde_Dparam: analytical derivative of matrix scaling var/cov matrix of shocks and measurement errors in transition equation with perturbation parameter w.r.t. param_identif
%       Detatilde_Dparam: analytical derivative of matrix scaling var/cov matrix of shocks and measurement errors in transition equation w.r.t. param_identif
%       Ddfstudt_Dparam: analytical derivative of degrees of freedom for t-distribution
%
% Outputs: m-file "filename_anal_deriv.m" 
%
% Based upon Schmitt-Grohé & Uribe (2013): Implementing Iskrev's Identifiabilty Test (http://www.columbia.edu/~mu2166/iskrev/iskrev.html)
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [] = anal_deriv_param_print2f(model,spec,approx,DSS_Dparam,Df_Donlyparam,Df_Dvars,Dgra_Donlyparam,Dgra_Dvars,Dhes_Donlyparam,Dhes_Dvars,DSigma_Dparam,Dsig_Dparam,Dsigetatilde_Dparam,Detatilde_Dparam,Ddfstudt_Dparam)

    %% Open m-file
    fid=fopen(['./models/',model,'/', model,'_spec',num2str(spec),'_approx',num2str(approx),'_anal_deriv.m'],'w');    
    
    %% PRINT JACOBIAN OF DSS_Dparam
        izero=(DSS_Dparam~=0);
        izeroi=find(izero);
        S=['nDSS_Dparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
        fprintf(fid,S);
        for i=1:length(izeroi);
         S2 = char(DSS_Dparam(izeroi(i)));
         S = [' nDSS_Dparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
         fprintf(fid,S);
        end
        S = ['nDSS_Dparam=reshape(nDSS_Dparam,[', num2str(size(izero)),']);\n'];  
        fprintf(fid,S); 
    
    %% PRINT JACOBIAN OF Dvars_Dparam = - inv(Df_Dvars)*Df_Donlyparam BY TERMS
        % Term 1: Derivative of f w.r.t param_identif only (without changes in steady-state)
        izero=(Df_Donlyparam~=0);
        izeroi=find(izero);
        S=['nDf_Donlyparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
        fprintf(fid,S);
        for i=1:length(izeroi);
         S2 = char(Df_Donlyparam(izeroi(i)));
         S = [' nDf_Donlyparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
         fprintf(fid,S);
        end
        S = ['nDf_Donlyparam=reshape(nDf_Donlyparam,[', num2str(size(izero)),']);\n'];  
        fprintf(fid,S);

        % Term 2: Derivative of f w.r.t steady-state variables (without changes in parameters)
         izero=(Df_Dvars~=0);
        izeroi=find(izero);
        S=['nDf_Dvars=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
        fprintf(fid,S);
        for i=1:length(izeroi);
         S2 = char(Df_Dvars(izeroi(i)));
         S = [' nDf_Dvars(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
         fprintf(fid,S);
        end
        S = ['nDf_Dvars=reshape(nDf_Dvars,[', num2str(size(izero)),']);\n'];  
        fprintf(fid,S);

        % Combine Term 1 and Term 2 to get Dvars_Dparam = - inv(Df_Dvars)*Df_Donlyparam
        S = ' nDvars_Dparam = - nDf_Dvars\\nDf_Donlyparam; \n';
        fprintf(fid,S);
    %%%
    
    %% PRINT JACOBIAN OF Dgra_Dparam BY TERMS
        % Term 1, i.e. derivative of gra w.r.t. param only (without changes in steady state)
        izero=(Dgra_Donlyparam~=0);
        izeroi=find(izero);
        S=['nDgra_Donlyparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
        fprintf(fid,S);
        for i=1:length(izeroi);
         S2 = char(Dgra_Donlyparam(izeroi(i)));
         S = [' nDgra_Donlyparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
         fprintf(fid,S);
        end
        S = ['nDgra_Donlyparam=reshape(nDgra_Donlyparam,[', num2str(size(izero)),']);\n'];  
        fprintf(fid,S);
        
        % Term 2, i.e. derivative of gra w.r.t. steady state (without changes in parameters)
        izero=(Dgra_Dvars~=0);
        izeroi=find(izero);
        S=['nDgra_Dvars=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
        fprintf(fid,S);
        for i=1:length(izeroi);
         S2 = char(Dgra_Dvars(izeroi(i)));
         S = [' nDgra_Dvars(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
         fprintf(fid,S);
        end
        S = ['nDgra_Dvars=reshape(nDgra_Dvars,[', num2str(size(izero)),']);\n'];  
        fprintf(fid,S);
        
        % Combine terms to get Dgra_Dparam
        S = ' nDgra_Dparam = nDgra_Donlyparam + nDgra_Dvars  * nDvars_Dparam; \n';
        fprintf(fid,S);
        S = ' clear nDgra_Donlyparam nDgra_Dvars; \n';
        fprintf(fid,S);        
    %%%

    %% PRINT JACOBIAN OF Dhes_Dparam BY TERMS
        % Term 1, i.e. derivative of hes w.r.t. param only (without changes in steady-state)
        izero=(Dhes_Donlyparam~=0);
        izeroi=find(izero);
        S=['nDhes_Donlyparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
        fprintf(fid,S);
        for i=1:length(izeroi);
         S2 = char(Dhes_Donlyparam(izeroi(i)));
         S = [' nDhes_Donlyparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
         fprintf(fid,S);
        end
        S = ['nDhes_Donlyparam=reshape(nDhes_Donlyparam,[', num2str(size(izero)),']);\n'];  
        fprintf(fid,S);
        
        % Term 2, i.e. derivative of hes w.r.t. steady-state (without changes in param)
        izero=(Dhes_Dvars~=0);
        izeroi=find(izero);
        S=['nDhes_Dvars=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
        fprintf(fid,S);
        for i=1:length(izeroi);
         S2 = char(Dhes_Dvars(izeroi(i)));
         S = [' nDhes_Dvars(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
         fprintf(fid,S);
        end
        S = ['nDhes_Dvars=reshape(nDhes_Dvars,[', num2str(size(izero)),']);\n'];  
        fprintf(fid,S); 

        % Combine terms to get Dhes_Dparam
        S = ' nDhes_Dparam = nDhes_Donlyparam + nDhes_Dvars  * nDvars_Dparam; \n';
        fprintf(fid,S);
        S = ' clear nDhes_Donlyparam nDhes_Dvars; \n';
        fprintf(fid,S);
    %%%
   
    %% PRINT JACOBIAN OF DSigma_Dparam
    izero=(DSigma_Dparam~=0);
    izeroi=find(izero);
    S=['nDSigma_Dparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi);
     S2 = char(DSigma_Dparam(izeroi(i)));
     S = [' nDSigma_Dparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDSigma_Dparam=reshape(nDSigma_Dparam,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
    %%%
    
    %% PRINT JACOBIAN OF Dsig_Dparam
    izero=(Dsig_Dparam~=0);
    izeroi=find(izero);
    S=['nDsig_Dparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi);
     S2 = char(Dsig_Dparam(izeroi(i)));
     S = [' nDsig_Dparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDsig_Dparam=reshape(nDsig_Dparam,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
    %%%

    %% PRINT JACOBIAN OF Dsigetatilde_Dparam
    izero=(Dsigetatilde_Dparam~=0);
    izeroi=find(izero);
    S=['nDsigetatilde_Dparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi);
     S2 = char(Dsigetatilde_Dparam(izeroi(i)));
     S = [' nDsigetatilde_Dparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDsigetatilde_Dparam=reshape(nDsigetatilde_Dparam,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
    %%%

    %% PRINT JACOBIAN OF Detatilde_Dparam
    izero=(Detatilde_Dparam~=0);
    izeroi=find(izero);
    S=['nDetatilde_Dparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi);
     S2 = char(Detatilde_Dparam(izeroi(i)));
     S = [' nDetatilde_Dparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDetatilde_Dparam=reshape(nDetatilde_Dparam,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
    %%%

    %% PRINT JACOBIAN OF Ddfstudt_Dparam
    izero=(Ddfstudt_Dparam~=0);
    izeroi=find(izero);
    S=['nDdfstudt_Dparam=zeros(',num2str(numel(izero)),',' num2str(1) ');\n'];
    fprintf(fid,S);
    for i=1:length(izeroi);
     S2 = char(Ddfstudt_Dparam(izeroi(i)));
     S = [' nDdfstudt_Dparam(' num2str(izeroi(i)) ',' num2str(1) ')= ' S2(1:end) ';  \n'];   
     fprintf(fid,S);
    end
    S = ['nDdfstudt_Dparam=reshape(nDdfstudt_Dparam,[', num2str(size(izero)),']);\n'];  
    fprintf(fid,S); 
    %%%

fclose(fid);
