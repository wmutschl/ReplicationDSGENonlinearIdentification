% This function 
% (1) Computes unique elements of the product moments of xi_t symbolically
% (2) Evaluates the unique elements of the product moments of xi_t for a Gaussian or t-distribution and writes them to files
% (3) Computes analytical derivatives of unique elements of product moments of xi_t for a Gaussian or t-distribution and writes them to files
% (4) Evaluates numerically the unique elments of the product moments of xi_t and its analytical derivatives
%
% Inputs:
%       nu: number of shocks, nx: number of states
%       Sigma, DSigma_Dparam: Covariance matrix of shocks and its Jacobian
%       E_xf_xf, DE_xf_xf_Dparam: Covariance matrix of states and its Jacobian
%       dfstudt, Ddfstudt_Dparam: degrees of freedom and its Jacobian for t-distribution
%       DSGE_Model: structure containing information about DSGE model
%       Ident_Test: structure containing information about which product momentes to compute
%       Settings: structure containing Settings%
% Calls:
%       duplication,triplication,quadruplication: Computes duplication, triplication and quadruplication matrices
%       ProdMom_inov_print2f: saves expressions to script files
%       prodmom: computes the product moment of u_{i_1}^{nu_1}u_{i_2}^{nu_2}...u_{i_m}^{nu_m}, where u_{i_j} are elements from u ~ N(0_n,Sigma)
%       allVL1
% Outputs:
%       nM2min, nDM2min_Dparam: numerically evaluated unique elements of E(kron(xi_t,xi_t)) and its Jacobian
%       nM3min, nDM3min_Dparam: numerically evaluated unique elements of E(kron(xi_t,xi_t,xi_t)) and its Jacobian
%       nM4min, nDM4min_Dparam: numerically evaluated unique elements of E(kron(xi_t,xi_t,xi_t,xi_t)) and its Jacobian
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)


function [nM2min,nM3min,nM4min,nDM2min_Dparam,nDM3min_Dparam,nDM4min_Dparam] = ProdMom_inov(nu,nx,Sigma,DSigma_Dparam,E_xf_xf,DE_xf_xf_Dparam,dfstudt,Ddfstudt_Dparam,DSGE_Model,Ident_Test,Settings)
% Initialize 
nM2min=[];nM3min=[];nM4min=[];nDM2min_Dparam=[];nDM3min_Dparam=[];nDM4min_Dparam=[];

%% Compute product moments symbolically depending on speed option
if strcmp(Settings.speed,'No Speed')
    % Create symbolic variables for u, make covariances positive
    u = sym('u',[nu 1]);
    SIGU = sym('sigu_%d%d',[nu nu]);
    assume(SIGU, 'positive');
    vecSIGUmin = SIGU(find(tril(ones(nu,nu))));
    SIGU= reshape(duplication(nu)*vecSIGUmin,nu,nu);
    % Create symbolic variables for xf, make covariances positive
    xf = sym('xf',[nx 1]);
    SIGX = sym('sigxf_%d%d',[nx nx]);
    assume(SIGX, 'positive');
    vecSIGXmin = SIGX(find(tril(ones(nx,nx))));
    SIGX= reshape(duplication(nx)*vecSIGXmin,nx,nx);
    SIG = [SIGU zeros(nu,nx); zeros(nx,nu) SIGX];
    % Create vector of auxiliary unique symbolic variables
    auxpar = [vecSIGUmin;vecSIGXmin];
    if strcmp(DSGE_Model.symbolic.distribution,'Student-t')
        v = sym('v');
        assume(v,'positive');
        auxpar = [auxpar;v];
    end
    % Create minimal xi_t vector
    ximin = [u; unique(kron(u,u),'stable')-vecSIGUmin; kron(u,xf)];
    nximin = nu + nu*(nu+1)/2 + nx*nu;

    %% Second-order moments symbolically
    if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
        fprintf('Get symbolic second-order moments of innovations\n');
        nunique2 = nximin*(nximin+1)/2;
        % % Get all combinations
        combos2 = flipud(allVL1(nximin, 2));% All integer permutations with sum criteria == 2
        S2 = cell(1,nximin);
        for ix = 1:nximin
            % string that contains unique product moments of each element in xi_t
            S2{ix} = ['ximin(',num2str(ix),')^i(',num2str(ix),')'];
        end
        S2 = strjoin(S2,'*');
        z2=sym(zeros(nunique2,1));
        reverseStr = '';
        for j = 1:nunique2
            if rem(j,50)== 0
                msg = sprintf('    Processed %d/%d', j, nunique2); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                    
            elseif j==nunique2
                msg = sprintf('    Processed %d/%d\n', j, nunique2); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                    
            end
            i = combos2(j,:); % powers
            % evaluate string of product moments of each element in xi_t
            z2(j,:) = eval(S2);
        end
        [z2min,ia2,ic2] = unique(z2,'stable'); %get rid of further duplicate elements
    else
        z2 =[]; z2min=[]; ia2=[]; ic2=[];
    end

    %% Third-order moments
    if Ident_Test.cumulants.third
        fprintf('Get symbolic third-order moments of innovations\n');
        nunique3 = nximin*(nximin+1)*(nximin+2)/6;
        % % Get all combinations
        combos3 = flipud(allVL1(nximin, 3)); % All integer permutations with sum criteria == 3
        S3 = cell(1,nximin);
        for ix = 1:nximin
            % string that contains unique product moments of each element in xi_t
            S3{ix} = ['ximin(',num2str(ix),')^i(',num2str(ix),')'];
        end
        S3 = strjoin(S3,'*');
        z3=sym(zeros(nunique3,1));
        reverseStr = '';
        for j = 1:nunique3
            if rem(j,50)== 0
               msg = sprintf('    Processed %d/%d', j, nunique3); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                    
            elseif j==nunique3
               msg = sprintf('    Processed %d/%d\n', j, nunique3); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                    
            end
            i = combos3(j,:); %powers
            % evaluate string of product moments of each element in xi_t
            z3(j,:) = eval(S3);
        end
        [z3min,ia3,ic3] = unique(z3,'stable');%get rid of further duplicate elements
    else
        z3 =[]; z3min=[]; ia3=[]; ic3=[];
    end

    %% Fourth-order moments
    if Ident_Test.cumulants.fourth
        fprintf('Get symbolic fourth-order moments of innovations\n');
        nunique4 = nximin*(nximin+1)*(nximin+2)*(nximin+3)/24;
        % % Get all combinations 
        combos4 = flipud(allVL1(nximin, 4)); % All integer permutations with sum criteria == 4
        S4 = cell(1,nximin);            
        for ix = 1:nximin
            % string that contains unique product moments of each element in xi_t
            S4{ix} = ['ximin(',num2str(ix),')^i(',num2str(ix),')'];                
        end
        S4 = strjoin(S4,'*');            
        z4=sym(zeros(nunique4,1));
        reverseStr = '';
        for j = 1:nunique4
            if rem(j,50)== 0
                msg = sprintf('    Processed %d/%d', j, nunique4); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                    
            elseif j==nunique4
                msg = sprintf('    Processed %d/%d\n', j, nunique4); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                    
            end
            i = combos4(j,:); %powers
            % evaluate string of product moments of each element in xi_t
            z4(j,:) = eval(S4);
        end        
        [z4min,ia4,ic4] = unique(z4,'stable');%get rid of further duplicate elements
    else
        z4 =[]; z4min=[]; ia4=[]; ic4=[];
    end

    %% Evaluate product moments analytically given theoretical moments from multivariate normal or t-distribution
    orderset = [2*(Ident_Test.cumulants.second||Ident_Test.cumulants.fourth) 3*Ident_Test.cumulants.third 4*Ident_Test.cumulants.fourth]; % which orders
    orderset(orderset==0) = []; M2=[]; M3=[]; M4=[]; %initialize
    for iorder = orderset
        fprintf('Evaluate product moments of xi analytically, order %d \n',iorder)
        z = expand(eval(['z',num2str(iorder),'min'])); %expand the symbolic expressions to have sums of products
        m = sym(zeros(size(z,1),1)); % initialize symbolic output vector
        reverseStr = ''; 
        for i=1:size(z,1)
            if rem(i,50)== 0
                msg = sprintf('    Processed %d/%d', i, size(z,1)); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
            elseif i==size(z,1)
                msg = sprintf('    Processed %d/%d\n', i, size(z,1)); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            STR = char(z(i)); % get individual entry as string
            % find summands: how many, + or -
            Sumsplus = strfind(STR,'+'); Sumsminus = strfind(STR,'-'); 
            Sums_idx1 = sort([1 (Sumsplus+2) (Sumsminus+2)]);
            Sums_idx2 = sort([(Sumsplus-1) (Sumsminus-2) length(STR)]);
            numSummands = numel(Sumsplus) + numel(Sumsminus)+1;        
            m_sum = sym(0); %initialize
            for iSum = 1:numSummands      % look at each summand individually      
                indx_mom = zeros(1,nu+nx); % intialize index that determines powers of E(u_1^{i_1}... u_{n_u}^{i_{n_u}) E(xf_1^{i_1}... xf_{n_x}^{i_{n_x})
                indx_rest = zeros(nu*(nu+1)/2,1); %initialize index that determines constant due to vec(SIGu)
                Str = STR(Sums_idx1(iSum):Sums_idx2(iSum));                
                isig=1;
                % Find indx_mom for u, and indx_rest 
                for iu=1:nu
                    str=['u',num2str(iu), '^'];
                    check = strfind(Str,str);
                    if ~isempty(check)
                        indx_mom(iu) = str2double(Str(check + length(str)));                
                        Str = strrep(Str,['u',num2str(iu), '^',Str(check + length(str))],'');
                    else
                        str(end) = [];
                        check = strfind(Str,str);
                        if ~isempty(check)
                            indx_mom(iu) = 1;
                            Str = strrep(Str,['u',num2str(iu)],'');
                        end        
                    end
                    %SIGMAS                
                    for ju=iu:nu                
                        str=['sigu_',num2str(ju),num2str(iu) '^'];
                        check = strfind(Str,str);
                        if ~isempty(check)
                            indx_rest(isig) = str2double(Str(check + length(str)));                
                            Str = strrep(Str,['sigu_',num2str(ju),num2str(iu), '^',Str(check + length(str))],'');
                        else
                            str(end) = [];
                            check = strfind(Str,str);
                            if ~isempty(check)
                                indx_rest(isig) = 1;
                                Str = strrep(Str,['sigu_',num2str(ju),num2str(iu),],'');
                            end        
                        end
                        isig=isig+1;
                    end
                end
                % Find indx_mom for xf
                for ix=1:nx
                    str=['xf',num2str(ix), '^'];
                    check = strfind(Str,str);
                    if ~isempty(check)
                        indx_mom(nu+ix) = str2double(Str(check + length(str)));
                        Str = strrep(Str,['xf',num2str(ix), '^',Str(check + length(str))],'');
                    else
                        str(end) = [];
                        check = strfind(Str,str);
                        if ~isempty(check)
                            indx_mom(nu+ix) = 1;
                            Str = strrep(Str,['xf',num2str(ix)],'');
                        end        
                    end        
                end
                
                % Find other constants, ie numbers
                Str = strrep(Str,'*','');
                const = strrep(Str,' ','');
                if isempty(str2double(const)) || isnan(str2double(const))
                    const = 1;
                else
                    const = str2double(const);
                end
                                
                if sum(indx_rest)==0
                    rest=1; % no terms belonging to powers of vec(SIGu)
                else
                    id_r = find(indx_rest~=0);
                    rest=1;
                    for r=1:length(id_r)
                        if strcmp(DSGE_Model.symbolic.distribution,'Student-t')                                
                            varstudt = v/(v-2);
                        else
                            varstudt = 1;
                        end
                        rest = rest*(varstudt*vecSIGUmin(id_r(r)))^indx_rest(id_r(r));
                    end
                end

                if sum(indx_mom)==0
                    mom=1;
                else                        
                    mom = prodmom(SIG,1:(nu+nx),indx_mom);
                    if strcmp(DSGE_Model.symbolic.distribution,'Student-t') && ~rem(sum(indx_mom),2)
                        dfu = mom_invGamma(v,sum(indx_mom(1:nu))/2);
                        dfx = mom_invGamma(v,sum(indx_mom(nu+1:end))/2);
                        mom = dfu*dfx*mom;
                    end
                end

                if iSum == 1 && strcmp(STR(1),'-')
                    sumsign = -1;
                elseif iSum > 1 && strcmp(STR(Sums_idx1(iSum)-2),'-')
                    sumsign = -1;
                else
                    sumsign = 1;
                end

                m_sum = m_sum + sumsign*const*rest*mom;        
            end % for end                
            m(i) = simplify(m_sum);
        end            
        eval(['M',num2str(iorder),'=m;'])
    end



    %% Analytical derivatives
    if strcmp(Settings.Derivative.type,'Analytical')
        if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
            fprintf('Computing analytical derivatives of product moments, order 2 \n');
            DM2_Dauxpar = simplify(jacobian(M2,auxpar));
        else
            DM2_Dauxpar = [];
        end
        if Ident_Test.cumulants.third
            fprintf('Computing analytical derivatives of product moments, order 3 \n');
            DM3_Dauxpar = simplify(jacobian(M3,auxpar));
        else
            DM3_Dauxpar = [];
        end
        if Ident_Test.cumulants.fourth
            fprintf('Computing analytical derivatives of product moments, order 4 \n');
            DM4_Dauxpar = simplify(jacobian(M4,auxpar));
        else
            DM4_Dauxpar = [];
        end
    else
        DM2_Dauxpar = []; 
        DM3_Dauxpar =[]; 
        DM4_Dauxpar=[];
    end

    %% Save to file
    ProdMom_inov_print2f(auxpar,DSGE_Model,Settings.approx,Settings.Derivative.type,Ident_Test,M2,M3,M4,DM2_Dauxpar,DM3_Dauxpar,DM4_Dauxpar,ic2,ic3,ic4);
end

%% Evaluate nummerically
% Set numerical values for auxiliary parameters in SIGU and SIGX and v
counti = 1; argu=zeros(nu*(nu+1)/2,1);
for j=1:nu
    for i=j:nu
        if strcmp(DSGE_Model.symbolic.distribution,'Student-t')
            argu(counti) = (dfstudt-2)/dfstudt*Sigma(i,j);
        else
            argu(counti) = Sigma(i,j);       
        end
        counti=counti+1;
    end
end
counti = 1; argx=zeros(nx*(nx+1)/2,1);
for j=1:nx
    for i=j:nx        
        if strcmp(DSGE_Model.symbolic.distribution,'Student-t')
            argx(counti) = (dfstudt-2)/dfstudt*E_xf_xf(i,j);
        else        
            argx(counti) = E_xf_xf(i,j);
        end
        counti=counti+1;
    end
end
arg = [argu;argx];
if strcmp(DSGE_Model.symbolic.distribution,'Student-t')
    arg = [arg;dfstudt];
end

if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth    
    fprintf('Evaluate script file prodmom2_num_eval from disk\n');
    try
       [nM2,ic2] = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom2_num_eval(arg)']);
    catch
        fprintf('   Script file is not yet on disk (strange bug), lets wait for 10 seconds otherwise just run again.\n');
        pause(10)
        fclose('all');
        [nM2,ic2] = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom2_num_eval(arg)']);
    end
    fprintf('Compute second-order cumulant of xi numerically \n');
    nM2min = sparse(nM2(ic2));
else
    nM2min = [];
end

if Ident_Test.cumulants.third
    fprintf('Evaluate script file prodmom3_num_eval from disk\n');
    try
       [nM3,ic3] = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom3_num_eval(arg)']);
    catch
        fprintf('   Script file is not yet on disk (strange bug), lets wait for 10 seconds otherwise just run again.\n');
        pause(10)
        fclose('all');
        [nM3,ic3] = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom3_num_eval(arg)']);
    end
    fprintf('Compute third-order cumulant of xi numerically \n');
    nM3min = sparse(nM3(ic3));    
else
    nM3min = [];
end

if Ident_Test.cumulants.fourth
    fprintf('Evaluate script file prodmom4_num_eval from disk\n');
    try
        [nM4,ic4] = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom4_num_eval(arg)']);        
    catch
        fprintf('   Script file is not yet on disk (strange bug), lets wait for 10 seconds otherwise just run again.\n');
        pause(10)
        fclose('all');
        [nM4,ic4] = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom4_num_eval(arg)']);
    end
    fprintf('Compute fourth-order cumulant of xi numerically \n');
    nM4min = sparse(nM4(ic4));
else
    nM4min = [];
end
  
if strcmp(Settings.Derivative.type,'Analytical')
    index= find(tril(ones(size(Sigma,1))));             
    nDvecSIGU_Dparam=DSigma_Dparam(index,:);
    index= find(tril(ones(size(E_xf_xf,1))));             
    nDvecSIGX_Dparam=DE_xf_xf_Dparam(index,:);    
    nDauxpar = sparse([nDvecSIGU_Dparam;nDvecSIGX_Dparam]);
    if strcmp(DSGE_Model.symbolic.distribution,'Student-t')
        nDauxpar = [nDauxpar;Ddfstudt_Dparam];
    end

    if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
        fprintf('Evaluate script file prodmom2_deriv from disk\n');
        try
            nDM2_Dauxpar = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom2_deriv(arg)']);            
        catch
            fprintf('   Script file is not yet on disk (strange bug), lets wait for 10 seconds otherwise just run again.\n');
            pause(10);
            fclose('all');
            nDM2_Dauxpar = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom2_deriv(arg)']);
        end
        fprintf('Compute analytical derivative of second-order cumulant of xi numerically \n');
        nDM2min_Dparam = sparse(nDM2_Dauxpar(ic2,:))*nDauxpar;
    end
    if Ident_Test.cumulants.third
        fprintf('Evaluate script file prodmom3_deriv from disk\n');
        try
            nDM3_Dauxpar = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom3_deriv(arg)']);
        catch
            fprintf('   Script file is not yet on disk (strange bug), lets wait for 10 seconds otherwise just run again.\n');
            pause(10);
            fclose('all');
            nDM3_Dauxpar = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom3_deriv(arg)']);
        end
        fprintf('Compute analytical derivative of third-order cumulant of xi numerically \n');
        nDM3min_Dparam = sparse(nDM3_Dauxpar(ic3,:))*nDauxpar;
    end
    if Ident_Test.cumulants.fourth
        fprintf('Evaluate script file prodmom4_deriv from disk\n');        
        try
            nDM4_Dauxpar = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom4_deriv(arg)']);
        catch
            fprintf('   Script file is not yet on disk (strange bug), lets wait for 10 seconds otherwise just run again.\n');
            pause(10);
            fclose('all');
            nDM4_Dauxpar = eval([DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx', num2str(Settings.approx), '_prodmom4_deriv(arg)']);
        end
        fprintf('Compute analytical derivative of fourth-order cumulant of xi numerically \n');
        nDM4min_Dparam = sparse(nDM4_Dauxpar(ic4,:))*nDauxpar;
    end
else
    nDM2min_Dparam=[]; nDM3min_Dparam=[]; nDM4min_Dparam=[];
end

function y = mom_invGamma(v,k)
    if k > 0
        a =v/2; b=v/2;
        divisor = a-1;
        for ii=2:k
            divisor = divisor*(a-ii);
        end
        y = b^k/divisor;
    elseif k== 0
        y=1;
    else
        error('Something went wrong with the invGamma moments')
    end    
end %mom_inv_gamma end

end % main function end
