% Computes numerical derivatives (up to second order) using two-sided finite difference method (also known as central differences)
% of expectation vector, solution matrices, var/cov matrix of shocks and measurement errors,
% of autocovariogram (for Iskrev's test), and of spectrum times its
% transpose (for Qu/Tkachenko's test)
%
% Inputs:
%       param_estim: local point at which to calculate derivatives
%       DSGE_Model: structure containing DSGE model independent of names
%       Settings: Structure containing settings about numerical derivatives and approximation order
%       Ident_Test: structure containing options for tests
%
% Calls: Embedded function GetObjects that calculates desired output objects at a local point xx
%
% Outputs:
%       Solut: Structure containing sparse numerical evaluated solution matrices, expectation, autocovariogram, and spectrum 
%       Deriv: Structure containing sparse numerical evaluated analytical derivative of solution matrices, expectation, autocovariogram, and spectrum w.r.t. param_identif (up to order 2)
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

function [Solut,Deriv] = NumericalDerivatives(param_estim,DSGE_Model,Settings,Ident_Test)
approx = Settings.approx;
nparam = length(DSGE_Model.param.identif);
DSGE_Model.numbers.nparam = nparam;
nx = DSGE_Model.numbers.nx;
nu = DSGE_Model.numbers.nu;
nparam_estim = length(DSGE_Model.param.estim);
param_names = DSGE_Model.param.names;

[Solut,Deriv] = EvaluateSparse(param_estim,DSGE_Model,approx,'Numerical');
%% Expectation Ed and Sigmax = reshape(E_xf_xf)
[Solut.Ed,~,Solut.E_xf_xf,~] = DerivExpectation(Solut,Deriv,DSGE_Model.numbers,'Numerical',approx);

%% Innovations (E_t(xi_t xi_t'))
%[Solut.Var_inov,~] = VarInnov(Solut.Sigma,[],Solut.E_xf_xf,[],Solut.vectorMom3,[],Solut.vectorMom4,[],nx,nu,'Numerical');
% Cumulants of Innovations xi_t=[u;kron(u,u)-vec(SIGU);kron(u,xf);kron(xf,u)]
[M2min,M3min,M4min,~,~,~] = ProdMom_inov(nu,nx,Solut.Sigma,[],Solut.E_xf_xf,[],Solut.dfstudt,[],DSGE_Model,Ident_Test,Settings);

Fxi = [speye(nu) spalloc(nu,nu*(nu+1)/2 + nu*nx,0);...
       spalloc(nu^2,nu,0) sparse(duplication(nu)) spalloc(nu^2,nu*nx,0);...
      spalloc(nu*nx,nu+nu*(nu+1)/2,0) speye(nu*nx); 
      spalloc(nx*nu,nu+nu*(nu+1)/2,0) sparse(commutation(nx,nu))];

% Save or load duplication matrix from file depending on speed setting
filename = ['./models/', DSGE_Model.shortname,'/',DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx',];
if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
    if strcmp(Settings.speed,'No Speed')
        fprintf('Compute duplication matrix\n');
        DPxi = sparse(duplication(DSGE_Model.numbers.nu+DSGE_Model.numbers.nu*(DSGE_Model.numbers.nu+1)/2+DSGE_Model.numbers.nu*DSGE_Model.numbers.nx));        
        try
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'DPxi','-append');
        catch
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'DPxi');
        end        
    else
        fprintf('Load duplication matrix for second-order cumulant\n');
        load([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'DPxi');
    end
else
    DPxi=[];
end

% Compute or load triplication matrix from file depending on speed setting
if Ident_Test.cumulants.third      
    if strcmp(Settings.speed,'No Speed')
        fprintf('Compute triplication matrix\n');
        TPxi = sparse(triplication(DSGE_Model.numbers.nu+DSGE_Model.numbers.nu*(DSGE_Model.numbers.nu+1)/2+DSGE_Model.numbers.nu*DSGE_Model.numbers.nx));
        try
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'TPxi','-append');
        catch
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'TPxi');
        end
    else
        fprintf('Load triplication matrix for third-order cumulant\n');
        load([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'TPxi');
    end
else
    TPxi = [];
end

% Save or load quadruplication matrix from file depending on speed setting
if Ident_Test.cumulants.fourth     
    if strcmp(Settings.speed,'No Speed')
        fprintf('Compute quadruplication matrix and its pseudoinverse\n');        
        [QPxi,QPxiinv] = quadruplication(DSGE_Model.numbers.nu+DSGE_Model.numbers.nu*(DSGE_Model.numbers.nu+1)/2+DSGE_Model.numbers.nu*DSGE_Model.numbers.nx);
        fprintf('Compute & save permutation matrix for fourth-order cumulant\n');
        % Permutation matrix for fourth-order cumulant
        p = DSGE_Model.numbers.nu+DSGE_Model.numbers.nu*(DSGE_Model.numbers.nu+1)/2+DSGE_Model.numbers.nu*DSGE_Model.numbers.nx;
        U = spalloc(p^3,p^3,p^3);
        for i=1:p^2
            for k=1:pp
                U((i-1)*p+k,(k-1)*p^2+i) = 1;        
            end
        end
        P = kron(speye(p),U);
        auxCUM4 = QPxiinv*(speye((DSGE_Model.numbers.nu+DSGE_Model.numbers.nu*(DSGE_Model.numbers.nu+1)/2+DSGE_Model.numbers.nu*DSGE_Model.numbers.nx)^4)+transpose(P)+P)*kron(DPxi,DPxi);
        try
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'auxCUM4','QPxi','-append');
        catch
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'auxCUM4','QPxi');
        end
    else
        fprintf('Load permutation matrix for fourth-order cumulant\n');
        load([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'auxCUM4','QPxi');
    end        
else
    auxCUM4 = []; QPxi=[];
end

%% Numerical derivatise options
% Set numerical differentiation step
tolnum = str2double(Settings.Derivative.options);

% Get selected parameters only
param_identif_numbers=[];
for p=1:nparam
    S = char(DSGE_Model.param.identif(p));
    param_identif_numbers = [param_identif_numbers str2double(S(6:end))];
end
% Numerical differentiation step for selected parameters only
for p=1:nparam_estim
    if isempty(find(param_identif_numbers==p))
        h(p)=0;
    else
        h(p) = tolnum*max(abs(param_estim(p)),1);            
    end
end
xh1=param_estim+h'; xh0=param_estim-h';
h=xh1-xh0; % h is now two times differentiation step

%% Computes two-sided finite difference Jacobians
for j=1:nparam;   
   xx = param_estim; % initialize local point
   % Solve model twice using pertubated local points xh1 and xh2
   xx(param_identif_numbers(j)) = xh1(param_identif_numbers(j));
   fprintf('...Calculating Numerical derivative for parameter %d out of %d left side \n',j,nparam)
   [Solut1,Iskrev_m1,QT_spec1] = GetObjects(xx,Fxi,DPxi,TPxi,QPxi,auxCUM4,DSGE_Model,Ident_Test,Settings,approx);
   fprintf('...Calculating Numerical derivative for parameter %d out of %d right side \n',j,nparam)
   xx(param_identif_numbers(j)) = xh0(param_identif_numbers(j)); 
   [Solut0,Iskrev_m0,QT_spec0] = GetObjects(xx,Fxi,DPxi,TPxi,QPxi,auxCUM4,DSGE_Model,Ident_Test,Settings,approx);
   % Approximate derivative using central differences
   DEd_Dparam(:,j) = (Solut1.Ed-Solut0.Ed)/h(param_identif_numbers(j));
   DSS_Dparam(:,j) = (Solut1.SS - Solut0.SS)/h(param_identif_numbers(j));
   Dprun_A_Dparam(:,j) = (Solut1.prun_A(:)-Solut0.prun_A(:))/h(param_identif_numbers(j));
   Dprun_B_Dparam(:,j) = (Solut1.prun_B(:)-Solut0.prun_B(:))/h(param_identif_numbers(j));
   Dprun_C_Dparam(:,j) = (Solut1.prun_C(:)-Solut0.prun_C(:))/h(param_identif_numbers(j));
   Dprun_D_Dparam(:,j) = (Solut1.prun_D(:)-Solut0.prun_D(:))/h(param_identif_numbers(j));
   Dprun_c_Dparam(:,j) = (Solut1.prun_c(:)-Solut0.prun_c(:))/h(param_identif_numbers(j));
   Dprun_d_Dparam(:,j) = (Solut1.prun_d(:)-Solut0.prun_d(:))/h(param_identif_numbers(j));
   DGAMMA2min_Dparam(:,j) = (Solut1.GAMMA2min(:)-Solut0.GAMMA2min(:))/h(param_identif_numbers(j));
   DGAMMA3min_Dparam(:,j) = (Solut1.GAMMA3min(:)-Solut0.GAMMA3min(:))/h(param_identif_numbers(j));
   DGAMMA4min_Dparam(:,j) = (Solut1.GAMMA4min(:)-Solut0.GAMMA4min(:))/h(param_identif_numbers(j));
   Iskrev_M(:,j) = (Iskrev_m1-Iskrev_m0)/h(param_identif_numbers(j));
   QT_dspec(:,(1+(j-1)*DSGE_Model.numbers.nd):j*DSGE_Model.numbers.nd) = (QT_spec1-QT_spec0)/h(param_identif_numbers(j));
    % Note: dspec is not the vectorized derivative (i.e. NOT Dspec_Dparam!)
    % It is a big matrix of dimension (N+1)*nd times nd*nparam           
    % For each of the (N+1) frequencies it consists of the non vectorized derivative 
    % of the nd times nd matrix Omega w.r.t. to the first parameter, then (next
    % to it) w.r.t. the second parameter and so on.
    % Then the same is done for the next frequency and stacked underneath
end

% Calculate objective function for Qu/Tkachenko's test
if (strcmp(Ident_Test.shortname,'QuTkachenko')||strcmp(Ident_Test.shortname,'rank'))
    if strcmp(Ident_Test.shortname,'QuTkachenko')
        N = str2double(Ident_Test.options{2}); % Number of subintervalls
    else
        N = str2double(Ident_Test.options{3}); % Number of subintervalls
    end
    QT_G=zeros(nparam,nparam); % Initialize G matrix
    j=1;
    nobs = DSGE_Model.numbers.nd;
    while j<=nparam; % For each parameter
        k=1;
        while k<=nparam; % For each parameter
           dOmegaj=QT_dspec(:,((j-1)*nobs+1):(j*nobs)); % dOmegaj is the nd times nd derivative of Omega w.r.t. thetaj, it is not vectorized
           dOmegak=QT_dspec(:,((k-1)*nobs+1):(k*nobs)); % dOmegak is the nd times nd derivative of Omega w.r.t. thetak, it is not vectorized
           f=1;
       while f<=(N+1); % For each frequency
    % Calculate the typical element of G
    QT_G(j,k)=QT_G(j,k)+trace(dOmegaj(((f-1)*nobs+1):(f*nobs),:)*dOmegak(((f-1)*nobs+1):(f*nobs),:));
    f=f+1; % Next frequency
       end
       k=k+1; % Next parameter
        end
        j=j+1; % Next parameter
    end
    QT_G=2*pi*QT_G./(N+1); % Approximate the integrand
else
    QT_G=0;
end

% Put everything into structure
Deriv.DEd_Dparam = DEd_Dparam;
Deriv.DSS_Dparam = DSS_Dparam;
Deriv.Dprun_A_Dparam = Dprun_A_Dparam;
Deriv.Dprun_B_Dparam = Dprun_B_Dparam;
Deriv.Dprun_C_Dparam = Dprun_C_Dparam;
Deriv.Dprun_D_Dparam = Dprun_D_Dparam;
Deriv.Dprun_c_Dparam = Dprun_c_Dparam;
Deriv.Dprun_d_Dparam = Dprun_d_Dparam;

Deriv.DGAMMA2min_Dparam = DGAMMA2min_Dparam;
Deriv.DGAMMA3min_Dparam = DGAMMA3min_Dparam;
Deriv.DGAMMA4min_Dparam = DGAMMA4min_Dparam;
Deriv.Iskrev_M = Iskrev_M;
Deriv.QT_G = QT_G;

%%
function [Sol,Iskrev_m,QT_S2] = GetObjects(xx,Fxi,DPxi,TPxi,QPxi,auxCUM4,DSGE_Model,Ident_Test,Settings,approx)
% Embedded function that calculates desired objects at a local point xx
%
% Inputs: 
%   xx: Local point for solving model
%   DSGE_Model: Structure that contains information about DSGE model independent of names
%   Ident_Test: Structure containing options about test
%   approx: order of approximation
%
% Calls:
%   EvaluateSparse: Evaluate model objects in sparse notation up to order 2
%   DerivExpectation: Calculate Expectation of observables up to order 2
%   ProdMom_inov: Cumulants of extended shock vector of innovations
    
% Outputs: 
%   Solut: Structure with solution matrices
%   Iskrev_m: Iskrev criteria
%   QT_spec: Qu and Tkachenko Criteria
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

%% Evaluate all model objects at steady-state numerically, solve model using sparse matrices
[Sol,~] = EvaluateSparse(xx,DSGE_Model,approx,'Numerical');
%% Expectation Ed and Sigmax = reshape(E_xf_xf)
[Sol.Ed,~,Sol.E_xf_xf,~] = DerivExpectation(Sol,[],DSGE_Model.numbers,'Numerical',approx);
%% Cumulants of Innovations xi_t=[u;kron(u,u)-vec(SIGU);kron(u,xf);kron(xf,u)]
[nM2min,nM3min,nM4min,~,~,~] = ProdMom_inov(DSGE_Model.numbers.nu,DSGE_Model.numbers.nx,Sol.Sigma,[],Sol.E_xf_xf,[],Sol.dfstudt,[],DSGE_Model,Ident_Test,Settings);
if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
    Sol.GAMMA2min = nM2min;
else
    Sol.GAMMA2min = [];
end
if Ident_Test.cumulants.third  
    Sol.GAMMA3min = nM3min;    
else
    Sol.GAMMA3min = [];
end
if Ident_Test.cumulants.fourth     
    Sol.GAMMA4min = nM4min -auxCUM4*kron(nM2min,nM2min);    
else
    Sol.GAMMA4min = [];
end

n_obs = DSGE_Model.numbers.nd;

if (strcmp(Ident_Test.shortname,'Iskrev')||strcmp(Ident_Test.shortname,'rank')) % Calculate cumulants of observables for Iskrev's test
    % Some precompuations
    T=str2double(Ident_Test.options{2});        % Get number of lags        
    BFxi = Sol.prun_B*Fxi; DFxi=Sol.prun_D*Fxi;
    if Ident_Test.cumulants.third || Ident_Test.cumulants.fourth
        AkronA = kron(Sol.prun_A,Sol.prun_A); CkronC = kron(Sol.prun_C,Sol.prun_C);
        BFxikronBFxi= kron(BFxi,BFxi); DFxikronDFxi= kron(DFxi,DFxi);        
        SelectMat2 = kron(Sol.SelectMat,Sol.SelectMat);
    end
    % Calculate zero-lag cumulant for second-order
    if Ident_Test.cumulants.second
        GAMMA2 = reshape(DPxi*Sol.GAMMA2min,size(Fxi,2),size(Fxi,2)); %Note GAMMA2 is a matrix, not vector
        BkronBGAMMA2 = BFxi*GAMMA2*transpose(BFxi);
        C2z0 = fastgensylv(-speye(size(Sol.prun_A,1)), Sol.prun_A, transpose(Sol.prun_A),BkronBGAMMA2,1e-17,1000);
        C2obs0 = Sol.SelectMat*Sol.prun_C*C2z0*transpose(Sol.SelectMat*Sol.prun_C) + Sol.SelectMat*DFxi*GAMMA2*transpose(Sol.SelectMat*DFxi);
        % Get index of same diagonal elements of zero-lag cumulant
        indexC2obs0= find(tril(ones(n_obs,n_obs)));
        C2obst = spalloc(n_obs^2,T-1,n_obs^2*(T-1)); 
    end
    % Calculate zero-lag cumulant for third-order
    if Ident_Test.cumulants.third        
        GAMMA3 = reshape(TPxi*Sol.GAMMA3min,size(Fxi,2)^2,size(Fxi,2)); %Note GAMMA3 is a matrix, not vector
        BkronBkronBGAMMA3 = BFxikronBFxi*GAMMA3*transpose(BFxi);
        C3z0 = fastgensylv(-speye(size(Sol.prun_A,1)^2), AkronA, transpose(Sol.prun_A),BkronBkronBGAMMA3,1e-17,1000);
        C3obs0 = SelectMat2*CkronC*C3z0*transpose(Sol.SelectMat*Sol.prun_C) + SelectMat2*DFxikronDFxi*GAMMA3*transpose(Sol.SelectMat*DFxi);
        % Get index of same diagonal elements of zero-lag cumulant
        indexC3obs0= find(repmat(tril(ones(n_obs,n_obs)),1,n_obs));
        C3obst = spalloc(n_obs^3,(T-1)*T/2,n_obs^3*((T-1)*T/2)); 
    end
    % Calculate zero-lag cumulant for fourth-order
    if Ident_Test.cumulants.fourth
        GAMMA4 = reshape(QPxi*Sol.GAMMA4min,size(Fxi,2)^2,size(Fxi,2)^2);%Note GAMMA4 is a matrix, not vector
        BkronBkronBkronBGAMMA4 = BFxikronBFxi*GAMMA4*transpose(BFxikronBFxi);
        C4z0 = fastgensylv(-speye(size(Sol.prun_A,1)^2), AkronA, transpose(AkronA),BkronBkronBkronBGAMMA4,1e-17,1000);
        C4obs0 = SelectMat2*CkronC*C4z0*transpose(SelectMat2*CkronC) + SelectMat2*DFxikronDFxi*GAMMA4*transpose(SelectMat2*DFxikronDFxi);
        % Get index of same diagonal elements of zero-lag cumulant
        indexC4obs0= find(repmat(tril(ones(n_obs,n_obs)),1,n_obs^2));    
        C4obst = spalloc(n_obs^4,(T-1)*T*(T+1)/6,n_obs^4*((T-1)*T*(T+1)/6));   
    end
    
    j1=1; j2=1; j3=1;
    At1=Sol.prun_A; 
    for t1=1:(T-1)
        if Ident_Test.cumulants.second
            C2zt = At1*C2z0;
            C2obst(:,j1) = vec(Sol.SelectMat*Sol.prun_C*C2zt*transpose(Sol.SelectMat*Sol.prun_C));
            j1=j1+1;
            At2=At1;
        else
            At2=At1;
        end
        for t2=t1:(T-1)
            if ~Ident_Test.cumulants.third && ~Ident_Test.cumulants.fourth
                break
            elseif ~Ident_Test.cumulants.third && Ident_Test.cumulants.fourth
                At3=At2;                      
            else
                At1kronAt2 = kron(At1,At2);
                C3zt = At1kronAt2*C3z0;
                C3obst(:,j2) = vec(SelectMat2*CkronC*C3zt*transpose(Sol.SelectMat*Sol.prun_C));
                j2=j2+1;
                At3=At2;
            end
                for t3=t2:(T-1)
                    if ~Ident_Test.cumulants.fourth
                        break
                    else
                        At2kronAt3 = kron(At2,At3);
                        IkronAt1T = kron(speye(size(At3,1)),transpose(At1));                        
                        C4zt = At2kronAt3*C4z0*IkronAt1T;
                        C4obst(:,j3) = vec(SelectMat2*CkronC*C4zt*transpose(SelectMat2*CkronC));
                        j3=j3+1;
                        At3=At3*Sol.prun_A;
                    end
                end
            At2=At2*Sol.prun_A;                
        end
        At1=At1*Sol.prun_A;
    end
    % Store cumulants of observables
    Iskrev_m = [];
    if Ident_Test.cumulants.second
        Iskrev_m=[Iskrev_m; vec(C2obs0(indexC2obs0)); C2obst(:)];
    end
    if Ident_Test.cumulants.third
        Iskrev_m=[Iskrev_m; vec(C3obs0(indexC3obs0)); C3obst(:)];
    end
    if Ident_Test.cumulants.fourth
        Iskrev_m=[Iskrev_m; vec(C4obs0(indexC4obs0)); C4obst(:)];
    end    
else
    Iskrev_m=[]; % If not Iskrev's test
end

if (strcmp(Ident_Test.shortname,'QuTkachenko')||strcmp(Ident_Test.shortname,'rank')) % Calculate derivative of spectra for Qu&Tkachenko's test    
    % Create vector of Fourier frequencies for approximation of the integral
    if strcmp(Ident_Test.shortname,'QuTkachenko')
        NN = str2double(Ident_Test.options{2}); % Number of subintervalls
    else
        NN = str2double(Ident_Test.options{3}); % Number of subintervalls
    end
    w=2*pi*(-(NN/2):1:(NN/2))'/NN; % Fourier frequencies    
    QT_S2=zeros(n_obs,n_obs,length(w)); % Initialize big matrix for power spectrum
    QT_S3=zeros(n_obs,n_obs,length(w)*(length(w)+1)/2); % Initialize big matrix for bispectrum
    QT_S4=zeros(n_obs,n_obs,length(w))*(length(w)+1)*(length(w)+2)/6; % Initialize big matrix for trispectrum
    for ii=1:length(w); %loop computes spectra
        z1=exp(-1i*w(ii)); %Lag Operator
        H1 = TransferFunction(z1,Sol.SelectMat,Sol.prun_A,Sol.prun_B,Sol.prun_C,Sol.prun_D); % z-transform
        QT_S2(:,:,ii)=(1/(2*pi))*H1*Fxi*GAMMA2*transpose(Fxi)*(H1'); % compute power spectrum
    end;
    QT_S2 = reshape(permute(QT_S2,[1 3 2]),DSGE_Model.numbers.nd*length(w),DSGE_Model.numbers.nd);
else
    QT_S2 = []; % If not Qu&Tkachenko's test
end

function H = TransferFunction(z,SelectMat,prun_A,prun_B,prun_C,prun_D)
    % Compute H         
    zIminusA =  (z*speye(size(prun_A,1)) - prun_A);
    zIminusAinv = zIminusA\speye(size(prun_A,1));
    H = SelectMat*(prun_D + prun_C*zIminusAinv*prun_B); % Transfer function    
end%transferfunction end

end% GetObjects end

end % Main function end