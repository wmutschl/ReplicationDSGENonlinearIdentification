% Computes analytical derivatives up to second order of expectation vector,
% solution matrices, cumulants and polyspectra
%
% Inputs:
%       param_estim: local point at which to calculate derivatives
%       DSGE_Model: structure containing DSGE model independent of names
%       Settings: structure containing settings
%       Ident_Test: structure containing settings for identification test
%
% Calls:
%       EvaluateSparse: Evaluate all model objects at steady-state numerically using sparse matrices
%       commutation: Calculates commutation matrix
%       DerivABCD:  Derivative of A*B*C*D w.r.t. param_identif
%       DerivExpectation: Derivative of expectation of observables
%       DerivXkronY: Derivative of kron(X,Y) w.r.t. param_identif
%       duplication,triplication,quadruplication: Computes auxiliary matrix for unique values of product moments of xi_t
%       ProdMom_inov: Compute analytically product moments and Jacobian of extended innovations vector xi_t
%       tracem: Matrix trace
%       TransferFunction: Embeded function that computes transfer function H(z)
%       vec: vectorize matrix
%
% Outputs:
%       Solut: Structure containing sparse numerical evaluated solution matrices, expectation, autocovariogram, and spectrum 
%       Deriv: Structure containing sparse numerical evaluated analytical derivative of solution matrices, expectation, autocovariogram, and spectrum w.r.t. param_identif (up to order 2)
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [Solut,Deriv] = AnalyticalDerivatives(param_estim,DSGE_Model,Settings,Ident_Test)
%% Evaluate all model objects at steady-state numerically using sparse matrices
[Solut,Deriv] = EvaluateSparse(param_estim,DSGE_Model,Settings.approx,'Analytical');
%% Some auxiliary precomputations
gra = Solut.gra; Dgra_Dparam= Deriv.Dgra_Dparam; 
hes = Solut.hes; Dhes_Dparam = Deriv.Dhes_Dparam; 
Sigma = Solut.Sigma; DSigma_Dparam = Deriv.DSigma_Dparam;
etatilde = Solut.etatilde; Detatilde_Dparam = Deriv.Detatilde_Dparam;
sig = Solut.sig; Dsig_Dparam = Deriv.Dsig_Dparam;
SelectMat = Solut.SelectMat; 
dfstudt = Solut.dfstudt; Ddfstudt_Dparam = Deriv.Ddfstudt_Dparam;
gv = Solut.gv; gvv = Solut.gvv; gSS = Solut.gSS;
hv = Solut.hv; hvv = Solut.hvv; hSS = Solut.hSS; 
prun_A = Solut.prun_A; prun_B = Solut.prun_B; prun_C = Solut.prun_C; prun_D = Solut.prun_D;
prun_c = Solut.prun_c; prun_d = Solut.prun_d;
SolM = Solut.solM; SolN = Solut.solN; SolQ = Solut.solQ; SolR = Solut.solR; 
SolS = Solut.solS; SolU = Solut.solU;

nv=DSGE_Model.numbers.nv; nx=DSGE_Model.numbers.nx; ny=DSGE_Model.numbers.ny; nu=DSGE_Model.numbers.nu;
nd=DSGE_Model.numbers.nd; nparam = size(Deriv.Dsig_Dparam,2);            
n=nv+ny;% Number of variables

%% Separate gradient
f1 = gra(:,1:nv);
Df1_Dparam= Dgra_Dparam(1:n*nv,:);
f2 = gra(:,nv+1:n);
Df2_Dparam= Dgra_Dparam((n*nv+1):(n^2),:);
f3 = gra(:,(n+1):(n+nv));
Df3_Dparam= Dgra_Dparam((n^2+1):(n^2+n*nv),:);
f4 = gra(:,(n+nv+1):end);
Df4_Dparam= Dgra_Dparam((n^2+n*nv+1):end,:);

%% Compute etaT eta_etaT and derivatives
etatildeT = transpose(etatilde);
DetatildeT_Dparam = sparse(commutation(size(etatilde))*Detatilde_Dparam);
etatilde_etatildeT = etatilde*transpose(etatilde);        
Detatilde_etatildeT_Dparam = DerivABCD(etatilde,Detatilde_Dparam,transpose(etatilde),DetatildeT_Dparam);
Solut.etatildeT = etatildeT;
Solut.etatilde_etatildeT = etatilde_etatildeT;
Deriv.DetatildeT_Dparam = DetatildeT_Dparam;
Deriv.Detatilde_etatildeT_Dparam = Detatilde_etatildeT_Dparam;

%% Construct Dhv_Dparam, DhvT_Dparam, Dgv_Dparam, DgvT_Dparam
F1=kron(transpose(hv),f2)+kron(speye(nv),f4);
F2=kron(speye(nv),f2*gv)+kron(speye(nv),f1);
F=-kron(transpose(hv)*transpose(gv),speye(n))*Df2_Dparam-kron(transpose(hv),speye(n))*Df1_Dparam -kron(transpose(gv),speye(n))*Df4_Dparam-Df3_Dparam;
Dgvhv_Dparam=[F1 F2]\F;
Dgv_Dparam = Dgvhv_Dparam(1:ny*nv,:);
Dhv_Dparam = Dgvhv_Dparam(ny*nv+1:end, :);
%Get transposes
DhvT_Dparam=sparse(commutation(size(hv)))*Dhv_Dparam;
DgvT_Dparam=sparse(commutation(size(gv)))*Dgv_Dparam;
Deriv.Dhv_Dparam = Dhv_Dparam; Deriv.DhvT_Dparam = DhvT_Dparam;
Deriv.Dgv_Dparam = Dgv_Dparam; Deriv.DgvT_Dparam = DgvT_Dparam;        

%% Construct Derivatives of second-order Solution matrices w.r.t param_identf
if Settings.approx==1 % Approximation to first order, set everything to zero
    Dgvv_Dparam=sparse(zeros(ny*nv^2,nparam)); 
    DgvvT_Dparam=sparse(zeros(ny*nv^2,nparam)); 
    DgSS_Dparam=sparse(zeros(ny,nparam)); 
    DgSST_Dparam=sparse(zeros(ny,nparam));     
    Dhvv_Dparam=sparse(zeros(nv^3,nparam)); 
    DhvvT_Dparam=sparse(zeros(nv^3,nparam)); 
    DhSS_Dparam=sparse(zeros(nv,nparam)); 
    DhSST_Dparam=sparse(zeros(nv,nparam));     
   
elseif Settings.approx==2 % Approximation to second order
    
    %% Construct Dgvv_Dparam, Dhvv_Dparam, DgvvT_Dparam, DhvvT_Dparam    
    
    % Derivative of Q1=kron(hv',f2,hv')+kron(eye(nv),f4,eye(nv))
    DSolQ1_Dparam = DerivXkronY(transpose(hv),DhvT_Dparam,kron(f2,transpose(hv)),DerivXkronY(f2,Df2_Dparam,transpose(hv),DhvT_Dparam))...
        + DerivXkronY(speye(nv),sparse(zeros(nv^2,nparam)),kron(f4,speye(nv)),DerivXkronY(f4,Df4_Dparam,speye(nv),sparse(zeros(nv^2,nparam))));
    % Derivative of Q2=kron(eye(nv),f1+f2*gv,eye(nv))
    Df2gv_Dparam = DerivABCD(f2,Df2_Dparam,gv,Dgv_Dparam);
    DSolQ2_Dparam = DerivXkronY(speye(nv),sparse(zeros(nv^2,nparam)),kron(f1+f2*gv,speye(nv)),DerivXkronY(f1+f2*gv,Df1_Dparam + Df2gv_Dparam,speye(nv),sparse(zeros(nv^2,nparam))));
    % Derivative of Q = [Q1 Q2]
    SolQinv = SolQ\speye(size(SolQ,1));
    DinvSolQ_Dparam=[]; % Using algorithm 1 of the paper       
    for i=1:nparam
        dSolQ = [reshape(DSolQ1_Dparam(:,i),n*nv^2,nv^2*ny) reshape(DSolQ2_Dparam(:,i),n*nv^2,nv^3)];
        dinvSolQ = -SolQinv*dSolQ*SolQinv;
        DinvSolQ_Dparam = [DinvSolQ_Dparam dinvSolQ(:)];
    end

    % Derivative of R=kron(eye(n),M')*hes*M with M=(hv, gv*hv,eye(nv),gv)'
    DSolM_Dparam =[];  % Using algorithm 1 of the paper    
    for i=1:nparam
        dSolM1 = reshape(Dhv_Dparam(:,i),size(hv));
        dSolM4 = reshape(Dgv_Dparam(:,i),size(gv));
        dSolM2 = gv*dSolM1 + dSolM4*hv;
        dSolM3 = zeros(nv);
        dSolM = [dSolM1; dSolM2; dSolM3; dSolM4];
        DSolM_Dparam = [DSolM_Dparam dSolM(:)];
    end
    DSolMT_Dparam = sparse(commutation(nv+ny+nv+ny,nv))*DSolM_Dparam;
    InKronSolMT = kron(speye(n),transpose(SolM));
    DInKronSolMT_Dparam = DerivXkronY(speye(n),sparse(zeros(n^2,nparam)),transpose(SolM),DSolMT_Dparam);
    DSolR_Dparam=DerivABCD(InKronSolMT,DInKronSolMT_Dparam,hes,Dhes_Dparam,SolM,DSolM_Dparam);

    % Derivative of [vec(gvv);vec(hvv)]=-Q^(-1)*vec(R)        
    Dgvvhvv_Dparam = (-1)*DerivABCD(SolQinv,DinvSolQ_Dparam,vec(SolR),DSolR_Dparam);
    Dgvv_Dparam = Dgvvhvv_Dparam(1:numel(gvv),:);
    DgvvT_Dparam = sparse(commutation(size(gvv)))*Dgvv_Dparam;
    Dhvv_Dparam = Dgvvhvv_Dparam(numel(gvv)+1:end,:);
    DhvvT_Dparam = sparse(commutation(size(hvv)))*Dhvv_Dparam;
    
    %% Construct DgSS_Dparam, DhSS_Dparam
    % Derivative of inv(S)=inv([S1 S2])=inv([f1+f2*gv f2+f4])
    SolSinv = SolS\speye(size(SolS,1));
    DSolS1_Dparam = Df1_Dparam + Df2gv_Dparam;
    DSolS2_Dparam = Df2_Dparam + Df4_Dparam;        
    DinvSolS_Dparam=[];  % Using algorithm 1 of the paper    
    for i=1:nparam
        dSolS = [reshape(DSolS1_Dparam(:,i),n,nv) reshape(DSolS2_Dparam(:,i),n,ny)];
        dinvSolS = -SolSinv*dSolS*SolSinv;
        DinvSolS_Dparam = [DinvSolS_Dparam dinvSolS(:)];
    end        

    % Derivative of U = f2*trm(U1) + trm(U2)
    % Derivative of U1 = kron(eye(ny),etatilde*etatilde')*gxx
    SolU1=kron(speye(ny),etatilde_etatildeT)*gvv;
    DSolU1_Dparam = DerivABCD(kron(speye(ny),etatilde_etatildeT),DerivXkronY(speye(ny),sparse(zeros(ny^2,nparam)),etatilde_etatildeT,Detatilde_etatildeT_Dparam),gvv,Dgvv_Dparam);
    % Derivative of U2 = kron(eye(n),N')*H*N*etatilde_etatildeT with N=(eye(nx),gx,zeros(n,nx))'        
    DSolN_Dparam = [];
    for i=1:nparam  % Using algorithm 1 of the paper    
        dSolN = [sparse(zeros(nv));reshape(Dgv_Dparam(:,i),size(gv));sparse(zeros(n,nv))];
        DSolN_Dparam=[DSolN_Dparam dSolN(:)];            
    end
    DSolNT_Dparam = sparse(commutation(2*n,nv))*DSolN_Dparam;     
    SolU2=kron(speye(n),transpose(SolN))*hes*SolN*etatilde_etatildeT;        
    DSolU2_Dparam = DerivABCD(kron(speye(n),transpose(SolN)),DerivXkronY(speye(n),sparse(zeros(n^2,nparam)),transpose(SolN),DSolNT_Dparam),hes,Dhes_Dparam,SolN,DSolN_Dparam,etatilde_etatildeT,Detatilde_etatildeT_Dparam);

    % Derivative of U = f2*trm(U1) + trm(U2)= trm(kron(eye(ny),etatilde*etatilde')*gxx) + trm(kron(eye(n),N')*H*N*eta_etaT with N=(eye(nx),gx,zeros(n,nx))'        )
    DSolU_Dparam=[];  % Using algorithm 1 of the paper    
    for i=1:nparam
        df2 = reshape(Df2_Dparam(:,i),size(f2));
        dtrmSolU1= sparse(tracem(reshape(DSolU1_Dparam(:,i),size(SolU1))));
        dtrmSolU2= sparse(tracem(reshape(DSolU2_Dparam(:,i),size(SolU2))));
        DSolU_Dparam = [DSolU_Dparam (df2*sparse(tracem(SolU1))+f2*dtrmSolU1+dtrmSolU2)];            
    end    

    % Derivative of [hSS;gSS]==-inv(S)*U        
    DhSS_gSS = -DerivABCD(SolSinv,DinvSolS_Dparam,SolU,DSolU_Dparam);
    DhSS_Dparam = DhSS_gSS(1:numel(hSS),:);
    DhSST_Dparam = sparse(commutation(size(hSS)))*DhSS_Dparam;
    DgSS_Dparam = DhSS_gSS(numel(hSS)+1:end,:);
    DgSST_Dparam = sparse(commutation(size(gSS)))*DgSS_Dparam;        
end

%% Store Derivatives of Solution matrices into structure
Deriv.Dhvv_Dparam = Dhvv_Dparam;
Deriv.DhvvT_Dparam = DhvvT_Dparam;
Deriv.Dgvv_Dparam = Dgvv_Dparam;
Deriv.DgvvT_Dparam = DgvvT_Dparam;
Deriv.DhSS_Dparam = DhSS_Dparam;
Deriv.DhSST_Dparam = DhSST_Dparam;
Deriv.DgSS_Dparam = DgSS_Dparam;
Deriv.DgSST_Dparam = DgSST_Dparam;

%% Construct Dprun_A_Dparam, Dprun_B_Dparam, Dprun_C_Dparam, Dprun_C_Dparam, Dprun_c_Dparam, Dprun_d_Dparam
hx = Solut.hv(Solut.ind.hx); hu = Solut.hv(Solut.ind.hu); hss=hSS(1:nx);
Hxx=Solut.hvv(Solut.ind.Hxx); Hxu=Solut.hvv(Solut.ind.Hxu); Hux=Solut.hvv(Solut.ind.Hux); Huu=Solut.hvv(Solut.ind.Huu);
Dhx_Dparam = Deriv.Dhv_Dparam(Solut.ind.hx,:); Dhu_Dparam = Deriv.Dhv_Dparam(Solut.ind.hu,:);
DHxx_Dparam=Deriv.Dhvv_Dparam(Solut.ind.Hxx,:); 
DHxu_Dparam=Deriv.Dhvv_Dparam(Solut.ind.Hxu,:); 
DHux_Dparam=Deriv.Dhvv_Dparam(Solut.ind.Hux,:); 
DHuu_Dparam=Deriv.Dhvv_Dparam(Solut.ind.Huu,:);
Dhss_Dparam=Deriv.DhSS_Dparam(1:nx,:);
Dhu_kron_hu = DerivXkronY(hu,Dhu_Dparam,hu,Dhu_Dparam);
Dhu_kron_hx = DerivXkronY(hu,Dhu_Dparam,hx,Dhx_Dparam);
Dhx_kron_hu = DerivXkronY(hx,Dhx_Dparam,hu,Dhu_Dparam);    
gx = Solut.gv(Solut.ind.gx); gu = Solut.gv(Solut.ind.gu);
Gxx=Solut.gvv(Solut.ind.Gxx); Gxu=Solut.gvv(Solut.ind.Gxu); Gux=Solut.gvv(Solut.ind.Gux); Guu=Solut.gvv(Solut.ind.Guu);
Dgx_Dparam = Deriv.Dgv_Dparam(Solut.ind.gx,:); Dgu_Dparam = Deriv.Dgv_Dparam(Solut.ind.gu,:);
DGxx_Dparam=Deriv.Dgvv_Dparam(Solut.ind.Gxx,:); 
DGxu_Dparam=Deriv.Dgvv_Dparam(Solut.ind.Gxu,:); 
DGux_Dparam=Deriv.Dgvv_Dparam(Solut.ind.Gux,:); 
DGuu_Dparam=Deriv.Dgvv_Dparam(Solut.ind.Guu,:);

Dprun_A_Dparam = []; Dprun_B_Dparam = []; Dprun_C_Dparam = []; Dprun_D_Dparam = []; Dprun_c_Dparam = []; Dprun_d_Dparam = [];
for i=1:nparam %Using algorithm 1 of the paper
    dprun_A = [reshape(Dhx_Dparam(:,i),size(hx)), sparse(zeros(nx,nx)),sparse(zeros(nx,nx^2));
               sparse(zeros(nx,nx)), reshape(Dhx_Dparam(:,i),size(hx)), 0.5*reshape(DHxx_Dparam(:,i),size(Hxx));
               sparse(zeros(nx*nx,nx)),sparse(zeros(nx*nx,nx)),reshape(DerivXkronY(hx,Dhx_Dparam(:,i),hx,Dhx_Dparam(:,i)),[nx^2,nx^2])];
    dprun_B = [reshape(Dhu_Dparam(:,i),size(hu)) sparse(zeros(nx,nu^2+nu*nx+nu*nx));...
              sparse(zeros(nx,nu)) 0.5*reshape(DHuu_Dparam(:,i),size(Huu)) 0.5*reshape(DHux_Dparam(:,i),size(Hux))  0.5*reshape(DHxu_Dparam(:,i),size(Hxu));...
              zeros(nx^2,nu) reshape(Dhu_kron_hu(:,i),[nx^2,nu^2]) reshape(Dhu_kron_hx(:,i),[nx^2,nu*nx]) reshape(Dhx_kron_hu(:,i),[nx^2,nx*nu])];
    dprun_C = [reshape(Dgx_Dparam(:,i),size(gx)), reshape(Dgx_Dparam(:,i),size(gx)),0.5*reshape(DGxx_Dparam(:,i),size(Gxx))];
    dprun_D = [reshape(Dgu_Dparam(:,i),size(gu)), 0.5*reshape(DGuu_Dparam(:,i),size(Guu)), 0.5*reshape(DGux_Dparam(:,i),size(Gux)) 0.5*reshape(DGxu_Dparam(:,i),size(Gxu))];
    dprun_c = [zeros(nx,1);
               reshape(0.5*(2*sig*hss*Dsig_Dparam(:,i) + sig^2*Dhss_Dparam(:,i)),size(hss)) + 0.5*(reshape(DHuu_Dparam(:,i),size(Huu))*vec(Sigma)+Huu*DSigma_Dparam(:,i));
               reshape(Dhu_kron_hu(:,i),[nx^2,nu^2])*vec(Sigma)+kron(hu,hu)*DSigma_Dparam(:,i)];
    dprun_d = reshape(0.5*(2*sig*gSS*Dsig_Dparam(:,i) + sig^2*DgSS_Dparam(:,i)),size(gSS))+0.5*(reshape(DGuu_Dparam(:,i),size(Guu))*vec(Sigma)+Guu*DSigma_Dparam(:,i));
    Dprun_A_Dparam = [Dprun_A_Dparam dprun_A(:)];
    Dprun_B_Dparam = [Dprun_B_Dparam dprun_B(:)];
    Dprun_C_Dparam = [Dprun_C_Dparam dprun_C(:)];
    Dprun_D_Dparam = [Dprun_D_Dparam dprun_D(:)];
    Dprun_c_Dparam = [Dprun_c_Dparam dprun_c(:)];
    Dprun_d_Dparam = [Dprun_d_Dparam dprun_d(:)];
end
% Store Derivatives of pruned-state-space into structure
Deriv.Dprun_A_Dparam = Dprun_A_Dparam;
Deriv.Dprun_B_Dparam = Dprun_B_Dparam;
Deriv.Dprun_C_Dparam = Dprun_C_Dparam;
Deriv.Dprun_D_Dparam = Dprun_D_Dparam; 
Deriv.Dprun_c_Dparam = Dprun_c_Dparam;
Deriv.Dprun_d_Dparam = Dprun_d_Dparam;

%% Analytical derivative of Expectation of observables (Ed)
[Ed,DEd_Dparam,E_xf_xf,DE_xf_xf_Dparam] = DerivExpectation(Solut,Deriv,DSGE_Model.numbers,'Analytical',Settings.approx);
Solut.Ed = Ed; Deriv.DEd_Dparam = DEd_Dparam;
Solut.E_xf_xf = E_xf_xf; Deriv.DE_xf_xf_Dparam = DE_xf_xf_Dparam;

%% Analytical Derivative of cumulants of Innovations xi_t=[u;kron(u,u)-vec(SIGU);kron(u,xf);kron(xf,u)]
[M2min,M3min,M4min,DM2min_Dparam,DM3min_Dparam,DM4min_Dparam] = ProdMom_inov(nu,nx,Sigma,DSigma_Dparam,E_xf_xf,DE_xf_xf_Dparam,dfstudt,Ddfstudt_Dparam,DSGE_Model,Ident_Test,Settings);

%% Save or load duplication matrix from file depending on speed setting
filename = ['./models/', DSGE_Model.shortname,'/',DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx',];
if Ident_Test.cumulants.second || Ident_Test.cumulants.fourth
    GAMMA2min = M2min; DGAMMA2min_Dparam = DM2min_Dparam;    
    if strcmp(Settings.speed,'No Speed')
        fprintf('Compute & save duplication matrix\n');
        DPxi = sparse(duplication(nu+nu*(nu+1)/2+nu*nx));
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
    GAMMA2min = []; DGAMMA2min_Dparam = [];
end
Deriv.DGAMMA2min_Dparam = DGAMMA2min_Dparam;

%% Save or load triplication matrix from file depending on speed setting
if Ident_Test.cumulants.third  
    GAMMA3min = M3min; DGAMMA3min_Dparam = DM3min_Dparam;
    Deriv.DGAMMA3min_Dparam = DGAMMA3min_Dparam;
    if strcmp(Settings.speed,'No Speed')
        fprintf('Compute & save triplication matrix\n');
        TPxi = sparse(triplication(nu+nu*(nu+1)/2+nu*nx));
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
    GAMMA3min = []; DGAMMA3min_Dparam = []; 
end
Deriv.DGAMMA3min_Dparam = DGAMMA3min_Dparam;

%% Save or load quadruplication matrix from file depending on speed setting
if Ident_Test.cumulants.fourth     
    if strcmp(Settings.speed,'No Speed')
        fprintf('Compute quadruplication matrix and its pseudoinverse\n');        
        [QPxi,QPxiinv] = quadruplication(nu+nu*(nu+1)/2+nu*nx);
        fprintf('Compute & save permutation matrix for fourth-order cumulant\n');
        % Permutation matrix for fourth-order cumulant
        p = nu+nu*(nu+1)/2+nu*nx;
        U = spalloc(p^3,p^3,p^3);
        for i=1:p^2
            for k=1:p
                U((i-1)*p+k,(k-1)*p^2+i) = 1;        
            end
        end
        P = kron(speye(p),U);
        auxCUM4 = QPxiinv*(speye((nu+nu*(nu+1)/2+nu*nx)^4)+transpose(P)+P)*kron(DPxi,DPxi);
        try
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'auxCUM4','QPxi','-append');
        catch
            save([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'auxCUM4','QPxi');
        end
    else
        fprintf('Load permutation matrix for fourth-order cumulant\n');
        load([filename,num2str(Settings.approx),'_prodmom_auxiliary'],'auxCUM4','QPxi');
    end    
    GAMMA4min = M4min -auxCUM4*kron(M2min,M2min);
    DGAMMA4min_Dparam = DM4min_Dparam - auxCUM4*DerivXkronY(M2min,DM2min_Dparam,M2min,DM2min_Dparam);
else
    GAMMA4min = []; DGAMMA4min_Dparam = [];
end
Deriv.DGAMMA4min_Dparam = DGAMMA4min_Dparam;

%% Iskrev's criteria
if (strcmp(Ident_Test.shortname,'Iskrev')||strcmp(Ident_Test.shortname,'rank'))
    % Construct cumulants analytically
    T=str2double(Ident_Test.options{2});        % Get number of lags    
    % Auxiliary matrix that selects unique elements in xi_t
    Fxi = [speye(nu) spalloc(nu,nu*(nu+1)/2 + nu*nx,0);...
           spalloc(nu^2,nu,0) sparse(duplication(nu)) spalloc(nu^2,nu*nx,0);...
          spalloc(nu*nx,nu+nu*(nu+1)/2,0) speye(nu*nx); 
          spalloc(nx*nu,nu+nu*(nu+1)/2,0) sparse(commutation(nx,nu))]; 
    
    % Calculate zero-lag cumulant for second-order
    if Ident_Test.cumulants.second                
        BFxi = prun_B*Fxi; DFxi=prun_D*Fxi;
        GAMMA2 = reshape(DPxi*GAMMA2min,size(Fxi,2),size(Fxi,2)); %Note GAMMA2 is a matrix, not vector
        BkronBGAMMA2 = BFxi*GAMMA2*transpose(BFxi);
        % Use fast General Sylvester algorithm 
        C2z0 = fastgensylv(-speye(size(prun_A,1)), prun_A, transpose(prun_A),BkronBGAMMA2, 1e-17,1000);
        for i=1:nparam
            dprun_A = reshape(Dprun_A_Dparam(:,i),size(prun_A));                        
            dprun_B = reshape(Dprun_B_Dparam(:,i),size(prun_B));
            dprun_C = reshape(Dprun_C_Dparam(:,i),size(prun_C));            
            dprun_D = reshape(Dprun_D_Dparam(:,i),size(prun_D));            
            dGAMMA2 = reshape(DPxi*DGAMMA2min_Dparam(:,i),size(Fxi,2),size(Fxi,2));
            dBkronBGAMMA2 = dprun_B*Fxi*GAMMA2*transpose(BFxi) + BFxi*dGAMMA2*transpose(BFxi) + BFxi*GAMMA2*transpose(dprun_B*Fxi);   
            dC2z0 = fastgensylv(-speye(size(prun_A,1)), prun_A, transpose(prun_A),prun_A*C2z0*transpose(dprun_A)+dprun_A*C2z0*transpose(prun_A)+dBkronBGAMMA2, 1e-17,1000);        
            dC2d0 = SelectMat*(dprun_C*C2z0*transpose(prun_C) + prun_C*dC2z0*transpose(prun_C) + prun_C*C2z0*transpose(dprun_C)...
                    + dprun_D*Fxi*GAMMA2*transpose(DFxi) + DFxi*dGAMMA2*transpose(DFxi) + DFxi*GAMMA2*transpose(dprun_D*Fxi))*transpose(SelectMat);
            DGAMMA2_Dparam(:,i) = dGAMMA2(:);
            DC2z0_Dparam(:,i) = dC2z0(:);
            DC2d0_Dparam(:,i) = dC2d0(:);                        
            clear dprun_A dprun_B dprun_C dprun_D dGAMMA2 dBkronBGAMMA2 dC2z0 dC2d0
        end
        % Get index of same diagonal elements of zero-lag cumulant
        indexC2d0= find(tril(ones(nd,nd)));        
        clear BkronBGAMMA2
    end
    
    % Calculate zero-lag cumulant for third-order  
    if Ident_Test.cumulants.third
        if ~Ident_Test.cumulants.second
            BFxi = prun_B*Fxi; DFxi=prun_D*Fxi;
        end
        SelectMat2 = kron(SelectMat,SelectMat);
        AkronA = kron(prun_A,prun_A);
        DAkronA_Dparam = DerivXkronY(prun_A,Dprun_A_Dparam,prun_A,Dprun_A_Dparam);
        BkronB = kron(prun_B,prun_B);
        DBkronB_Dparam = DerivXkronY(prun_B,Dprun_B_Dparam,prun_B,Dprun_B_Dparam);
        CkronC = kron(prun_C,prun_C);
        DCkronC_Dparam = DerivXkronY(prun_C,Dprun_C_Dparam,prun_C,Dprun_C_Dparam);
        DkronD = kron(prun_D,prun_D);
        DDkronD_Dparam = DerivXkronY(prun_D,Dprun_D_Dparam,prun_D,Dprun_D_Dparam);        
        FxikronFxi = kron(Fxi,Fxi);
        BFxikronBFxi= kron(BFxi,BFxi); DFxikronDFxi= kron(DFxi,DFxi);        
        GAMMA3 = reshape(TPxi*GAMMA3min,size(Fxi,2)^2,size(Fxi,2)); %Note GAMMA3 is a matrix, not vector
        BkronBkronBGAMMA3 = BFxikronBFxi*GAMMA3*transpose(BFxi);
        % Use fast General Sylvester algorithm 
        C3z0 = fastgensylv(-speye(size(prun_A,1)^2), AkronA, transpose(prun_A),BkronBkronBGAMMA3, 1e-17,1000);
        for i=1:nparam
            dprun_A = reshape(Dprun_A_Dparam(:,i),size(prun_A));                        
            dprun_B = reshape(Dprun_B_Dparam(:,i),size(prun_B));
            dprun_C = reshape(Dprun_C_Dparam(:,i),size(prun_C));            
            dprun_D = reshape(Dprun_D_Dparam(:,i),size(prun_D));
            dAkronA = reshape(DAkronA_Dparam(:,i),size(prun_A,1)^2,size(prun_A,2)^2);
            dBkronB = reshape(DBkronB_Dparam(:,i),size(prun_B,1)^2,size(prun_B,2)^2);
            dCkronC = reshape(DCkronC_Dparam(:,i),size(prun_C,1)^2,size(prun_C,2)^2);
            dDkronD = reshape(DDkronD_Dparam(:,i),size(prun_D,1)^2,size(prun_D,2)^2);
            dGAMMA3 = reshape(TPxi*DGAMMA3min_Dparam(:,i),size(Fxi,2)^2,size(Fxi,2));
            dBkronBkronBGAMMA3 = dBkronB*FxikronFxi*GAMMA3*transpose(BFxi) + BFxikronBFxi*dGAMMA3*transpose(BFxi) + BFxikronBFxi*GAMMA3*transpose(dprun_B*Fxi);   
            dC3z0 = fastgensylv(-speye(size(prun_A,1)^2), AkronA, transpose(prun_A),AkronA*C3z0*transpose(dprun_A)+dAkronA*C3z0*transpose(prun_A)+dBkronBkronBGAMMA3, 1e-17,1000);
            dC3d0 = SelectMat2*(dCkronC*C3z0*transpose(prun_C) + CkronC*dC3z0*transpose(prun_C) + CkronC*C3z0*transpose(dprun_C)...
                    + dDkronD*FxikronFxi*GAMMA3*transpose(DFxi) + DFxikronDFxi*dGAMMA3*transpose(DFxi) + DFxikronDFxi*GAMMA3*transpose(dprun_D*Fxi) )*transpose(SelectMat);
            DGAMMA3_Dparam(:,i) = dGAMMA3(:);
            DC3z0_Dparam(:,i) = dC3z0(:);
            DC3d0_Dparam(:,i) = dC3d0(:);
            clear dprun_A dprun_B dprun_C dprun_D dAkronA dBkronB dDkronD dGAMMA3 dBkronBkronBGAMMA3 dC3z0 dC3d0
        end
        % Get index of same diagonal elements of zero-lag cumulant
        indexC3d0= find(repmat(tril(ones(nd,nd)),1,nd));
        clear BkronBkronBGAMMA3 
    end
    
    % Calculate zero-lag cumulant for fourth-order
    if Ident_Test.cumulants.fourth
        if ~Ident_Test.cumulants.second
            BFxi = prun_B*Fxi; DFxi=prun_D*Fxi;
        end
        if ~Ident_Test.cumulants.third
            SelectMat2 = kron(SelectMat,SelectMat);
            AkronA = kron(prun_A,prun_A);
            DAkronA_Dparam = DerivXkronY(prun_A,Dprun_A_Dparam,prun_A,Dprun_A_Dparam);
            DBkronB_Dparam = DerivXkronY(prun_B,Dprun_B_Dparam,prun_B,Dprun_B_Dparam);
            CkronC = kron(prun_C,prun_C);
            DCkronC_Dparam = DerivXkronY(prun_C,Dprun_C_Dparam,prun_C,Dprun_C_Dparam);
            DDkronD_Dparam = DerivXkronY(prun_D,Dprun_D_Dparam,prun_D,Dprun_D_Dparam);        
            FxikronFxi = kron(Fxi,Fxi);
            BFxikronBFxi= kron(BFxi,BFxi); DFxikronDFxi= kron(DFxi,DFxi);
        end
        GAMMA4 = reshape(QPxi*GAMMA4min,size(Fxi,2)^2,size(Fxi,2)^2);%Note GAMMA4 is a matrix, not vector
        BkronBkronBkronBGAMMA4 = BFxikronBFxi*GAMMA4*transpose(BFxikronBFxi);
        % Use fast General Sylvester algorithm 
        C4z0 = fastgensylv(-speye(size(prun_A,1)^2), AkronA, transpose(AkronA),BkronBkronBkronBGAMMA4, 1e-17,1000);
        for i=1:nparam
            dAkronA = reshape(DAkronA_Dparam(:,i),size(prun_A,1)^2,size(prun_A,2)^2);
            dBkronB = reshape(DBkronB_Dparam(:,i),size(prun_B,1)^2,size(prun_B,2)^2);
            dCkronC = reshape(DCkronC_Dparam(:,i),size(prun_C,1)^2,size(prun_C,2)^2);
            dDkronD = reshape(DDkronD_Dparam(:,i),size(prun_D,1)^2,size(prun_D,2)^2);
            dGAMMA4 = reshape(QPxi*DGAMMA4min_Dparam(:,i),size(Fxi,2)^2,size(Fxi,2)^2);
            dBkronBkronBkronBGAMMA4 = dBkronB*FxikronFxi*GAMMA4*transpose(BFxikronBFxi) + BFxikronBFxi*dGAMMA4*transpose(BFxikronBFxi) + BFxikronBFxi*GAMMA4*transpose(dBkronB*FxikronFxi);               
            dC4z0 = fastgensylv(-speye(size(prun_A,1)^2), AkronA, transpose(AkronA),AkronA*C4z0*transpose(dAkronA)+dAkronA*C4z0*transpose(AkronA)+dBkronBkronBkronBGAMMA4, 1e-17,1000);
            dC4d0 = SelectMat2*(dCkronC*C4z0*transpose(CkronC) + CkronC*dC4z0*transpose(CkronC) + CkronC*C4z0*transpose(dCkronC)...
                    + dDkronD*FxikronFxi*GAMMA4*transpose(DFxikronDFxi) + DFxikronDFxi*dGAMMA4*transpose(DFxikronDFxi) + DFxikronDFxi*GAMMA4*transpose(dDkronD*FxikronFxi) )*transpose(SelectMat2);
            DGAMMA4_Dparam(:,i) = dGAMMA4(:);
            DC4z0_Dparam(:,i) = dC4z0(:);
            DC4d0_Dparam(:,i) = dC4d0(:);
            clear dAkronA dBkronB dDkronD dGAMMA4 dBkronBkronBkronBGAMMA4 dC4z0 dC4d0
        end
        % Get index of same diagonal elements of zero-lag cumulant
        indexC4d0= find(repmat(tril(ones(nd,nd)),1,nd^2));    
        clear BkronBkronBkronBGAMMA4
    end
    
    % Initialize storage for cumulants
    DC2dt_Dparam =sparse(zeros((T-1)*nd^2,nparam));
    DC3dt_Dparam =sparse(zeros((T-1)*T/2*nd^3,nparam)); 
    DC4dt_Dparam =sparse(zeros((T-1)*T*(T+1)/6*nd^4,nparam));    
    for j=1:nparam
        reverseStr = '';
        fprintf('Compute Cumulantogram for parameter j=%d of %d\n',j,nparam);
        dprun_A = reshape(Dprun_A_Dparam(:,j),size(prun_A));
        At1=prun_A; 
        dAt1=dprun_A;
        dprun_C = reshape(Dprun_C_Dparam(:,j),size(prun_C));
        if Ident_Test.cumulants.second
            dC2z0 = reshape(DC2z0_Dparam(:,j), size(C2z0));
        end
        if Ident_Test.cumulants.third
            dC3z0 = reshape(DC3z0_Dparam(:,j), size(C3z0));
        end
        if Ident_Test.cumulants.fourth            
            dC4z0 = reshape(DC4z0_Dparam(:,j), size(C4z0));
        end
        % Compute derivative of second-order Cumulants of observables
        i2=1; i3=1; i4=1;
        for t1=1:T-1
            if rem(t1,10) == 0
                msg = sprintf('   t1=%d of %d',t1,T-1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
            elseif t1==T-1
                msg = sprintf('   t1=%d of %d\n',t1,T-1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                
            end
            if Ident_Test.cumulants.second
                C2zt = At1*C2z0;
                dC2zt = dAt1*C2z0 + At1*dC2z0;
                dC2dt = SelectMat*(dprun_C*C2zt*transpose(prun_C) + prun_C*dC2zt*transpose(prun_C) + prun_C*C2zt*transpose(dprun_C))*transpose(SelectMat);
                DC2dt_Dparam(nd^2*(i2-1)+1:nd^2*i2,j) = dC2dt(:);
                i2=i2+1;
            end
            At2=At1; dAt2=dAt1;
            % Compute derivative of third-order Cumulants of observables
            for t2=t1:T-1
                if ~Ident_Test.cumulants.third && ~Ident_Test.cumulants.fourth
                    break
                elseif ~Ident_Test.cumulants.third && Ident_Test.cumulants.fourth
                    At3=At2; dAt3=dAt2;
                elseif Ident_Test.cumulants.third
                    if rem(t2,10) == 0
                        msg = sprintf('   t1=%d of %d, t2=%d of %d each, with t1<t2',t1,T-1,t2,T-1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    elseif t1==T-1 && t2==T-1
                        msg = sprintf('   t1=%d of %d, t2=%d of %d each, with t1<t2\n',t1,T-1,t2,T-1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                    At1kronAt2 = kron(At1,At2);
                    dAt1kronAt2 = kron(dAt1,At2) + kron(At1,dAt2);                    
                    C3zt = At1kronAt2*C3z0;
                    dC3zt = dAt1kronAt2*C3z0 + At1kronAt2*dC3z0;
                    dC3dt = SelectMat2*(dCkronC*C3zt*transpose(prun_C) + CkronC*dC3zt*transpose(prun_C) + CkronC*C3zt*transpose(dprun_C))*transpose(SelectMat);
                    DC3dt_Dparam(nd^3*(i3-1)+1:nd^3*i3,j) = dC3dt(:);                    
                    i3=i3+1;
                    At3=At2; dAt3=dAt2;                    
                end                
                % Compute derivative of fourth-order Cumulants of observables
                for t3=t2:T-1
                    if ~Ident_Test.cumulants.fourth
                        break
                    else
                        if rem(t3,10) == 0
                            msg = sprintf('   t1=%d of %d, t2=%d of %d, t3=%d of %d each, with t1<t2<t3',t1,T-1,t2,T-1,t3,T-1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));                        
                        elseif t1==T-1 && t2==T-1 && t3==T-1
                            msg = sprintf('   t1=%d of %d, t2=%d of %d, t3=%d of %d each, with t1<t2<t3\n',t1,T-1,t2,T-1,t3,T-1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                        end
                        At2kronAt3 = kron(At2,At3);
                        dAt2kronAt3 = kron(dAt2,At3) + kron(At2,dAt3);
                        IkronAt1T = kron(speye(size(At3,1)),transpose(At1));
                        dIkronAt1T = kron(speye(size(At3,1)),transpose(dAt1));
                        C4zt = At2kronAt3*C4z0*IkronAt1T;
                        dC4zt = dAt2kronAt3*C4z0*IkronAt1T + At2kronAt3*dC4z0*IkronAt1T + At2kronAt3*C4z0*dIkronAt1T;
                        dC4dt = SelectMat2*(dCkronC*C4zt*transpose(CkronC) + CkronC*dC4zt*transpose(CkronC) + CkronC*C4zt*transpose(dCkronC))*transpose(SelectMat2);
                        DC4dt_Dparam(nd^4*(i4-1)+1:nd^4*i4,j) = dC4dt(:);
                        % Calculate new DAt3_Dparam and At3
                        dAt3=At3*dprun_A + dAt3*prun_A;
                        At3=At3*prun_A;
                        i4=i4+1;                        
                    end
                end
                % Calculate new DAt2_Dparam and At2
                dAt2=At2*dprun_A + dAt2*prun_A;
                At2=At2*prun_A;                
            end
            % Calculate new DAt1_Dparam and At1
            dAt1=At1*dprun_A + dAt1*prun_A;
            At1=At1*prun_A;
        end        
    end
    
    % Store cumulants of observables
    Iskrev_M = [];
    if Ident_Test.cumulants.second
        Iskrev_M=[Iskrev_M; DC2d0_Dparam(indexC2d0,:); DC2dt_Dparam];
    end
    if Ident_Test.cumulants.third
        Iskrev_M=[Iskrev_M; DC3d0_Dparam(indexC3d0,:); DC3dt_Dparam];
    end
    if Ident_Test.cumulants.fourth
        Iskrev_M=[Iskrev_M; DC4d0_Dparam(indexC4d0,:); DC4dt_Dparam];
    end
else
    Iskrev_M=[]; % If not Iskrev's test
end
Deriv.Iskrev_M = Iskrev_M;

%% Qu and Tkachenko's criteria
if (strcmp(Ident_Test.shortname,'QuTkachenko')||strcmp(Ident_Test.shortname,'rank'))
    % Construct G analytically
    fprintf('Compute Polyspectra for frequencies\n')
    % Create vector of Fourier frequencies for approximation of the integral
    if strcmp(Ident_Test.shortname,'QuTkachenko')
        N=str2double(Ident_Test.options{2}); % Subintervalls
    else
        N=str2double(Ident_Test.options{3}); % Subintervalls
    end
    % Create Fourier frequencies
    w=2*pi*(-(N/2):1:(N/2))'/N;
    % Calculate zero-lag cumulants and its derivative
    if ~strcmp(Ident_Test.shortname,'rank')
        % Auxiliary matrix that selects unique elements in xi_t
        Fxi = [speye(nu) spalloc(nu,nu*(nu+1)/2 + nu*nx,0);...
                spalloc(nu^2,nu,0) sparse(duplication(nu)) spalloc(nu^2,nu*nx,0);...
                spalloc(nu*nx,nu+nu*(nu+1)/2,0) speye(nu*nx); 
                spalloc(nx*nu,nu+nu*(nu+1)/2,0) sparse(commutation(nx,nu))]; 

        if Ident_Test.cumulants.second            
            GAMMA2 = reshape(DPxi*GAMMA2min,size(Fxi,2),size(Fxi,2)); %Note GAMMA2 is a matrix, not a vector
        end
        if Ident_Test.cumulants.third
            SelectMat2 = kron(SelectMat,SelectMat);     
            GAMMA3 = reshape(TPxi*GAMMA3min,size(Fxi,2)^2,size(Fxi,2)); %Note GAMMA3 is a matrix, not a vector
        end
        if Ident_Test.cumulants.fourth
            GAMMA4 = reshape(QPxi*GAMMA4min,size(Fxi,2)^2,size(Fxi,2)^2); %Note GAMMA4 is a matrix, not a vector
        end

        for i=1:nparam  
            if Ident_Test.cumulants.second
                dGAMMA2 = reshape(DPxi*DGAMMA2min_Dparam(:,i),size(Fxi,2),size(Fxi,2));
                DGAMMA2_Dparam(:,i) = dGAMMA2(:);
                clear dGAMMA2
            end
            if Ident_Test.cumulants.third
                dGAMMA3 = reshape(TPxi*DGAMMA3min_Dparam(:,i),size(Fxi,2)^2,size(Fxi,2));
                DGAMMA3_Dparam(:,i) = dGAMMA3(:);
                clear dGAMMA3
            end
            if Ident_Test.cumulants.fourth
                dGAMMA4 = reshape(QPxi*DGAMMA4min_Dparam(:,i),size(Fxi,2)^2,size(Fxi,2)^2);
                DGAMMA4_Dparam(:,i) = dGAMMA4(:);
                clear dGAMMA4
            end             
        end
    end
    FxikronFxi = kron(Fxi,Fxi);
    if Ident_Test.cumulants.second
        GAMMA2full = Fxi*GAMMA2*transpose(Fxi);
    end
    if Ident_Test.cumulants.third
        GAMMA3full = FxikronFxi*GAMMA3*transpose(Fxi);
    end
    if Ident_Test.cumulants.fourth
        GAMMA4full = FxikronFxi*GAMMA4*transpose(FxikronFxi);
    end
    for j=1:nparam
        if Ident_Test.cumulants.second
        DGAMMA2full_Dparam(:,j) = vec(Fxi*reshape(DGAMMA2_Dparam(:,j),size(Fxi,2),size(Fxi,2))*transpose(Fxi));
        end
        if Ident_Test.cumulants.third
            DGAMMA3full_Dparam(:,j) = vec(FxikronFxi*reshape(DGAMMA3_Dparam(:,j),size(Fxi,2)^2,size(Fxi,2))*transpose(Fxi));
        end
        if Ident_Test.cumulants.fourth
            DGAMMA4full_Dparam(:,j) = vec(FxikronFxi*reshape(DGAMMA4_Dparam(:,j),size(Fxi,2)^2,size(Fxi,2)^2)*transpose(FxikronFxi));
        end
    end
    G2 = zeros(nparam,nparam); % Initialize objective function for power spectrum
    G3 = zeros(nparam,nparam); % Initialize objective function for bispectrum
    G4 = zeros(nparam,nparam); % Initialize objective function for trispectrum
    reverseStr = '';
    for l1=1:length(w); %loop computes analytical derivative of spectra
        z1=exp(-1i*w(l1)); % Use Fourier transform for lag operator
        [Hz1,DHz1_Dparam,DHz1cT_Dparam] = TransferFunction(z1,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam);
        % Compute derivative of power spectrum
        if Ident_Test.cumulants.second
            DOmega2_Dparam = (1/(2*pi))*DerivABCD(Hz1,DHz1_Dparam,GAMMA2full,DGAMMA2full_Dparam,Hz1',DHz1cT_Dparam);
            %G2 = G2 + transpose(commutation(nd,nd)*DOmega2_Dparam)*DOmega2_Dparam;
            G2 = G2 + DOmega2_Dparam'*DOmega2_Dparam;
        end
        for l2=l1:length(w); %loop computes analytical derivative of spectra and stacks them in a big matrix
            if ~Ident_Test.cumulants.third && ~Ident_Test.cumulants.fourth
                break
            elseif ~Ident_Test.cumulants.third && Ident_Test.cumulants.fourth
                z2=exp(-1i*w(l2)); % Use Fourier transform for lag operator
                [Hz2,DHz2_Dparam,DHz2cT_Dparam] = TransferFunction(z2,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam);
            elseif Ident_Test.cumulants.third
                z2=exp(-1i*w(l2)); % Use Fourier transform for lag operator
                [Hz2,DHz2_Dparam,DHz2cT_Dparam] = TransferFunction(z2,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam);
                Hz1kronHz2 = kron(Hz1,Hz2);
                DHz1kronHz2_Dparam = DerivXkronY(Hz1,DHz1_Dparam,Hz2,DHz2_Dparam);
                [Hz1z2,DHz1z2_Dparam,DHz1z2cT_Dparam] = TransferFunction(z1*z2,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam);
                % Compute derivative of bispectrum
                DOmega3_Dparam = (1/(4*pi^2))*DerivABCD(Hz1kronHz2,DHz1kronHz2_Dparam,GAMMA3full,DGAMMA3full_Dparam,Hz1z2',DHz1z2cT_Dparam);
                %G3 = G3 + transpose(commutation(nd^2,nd)*DOmega3_Dparam)*DOmega3_Dparam;                
                G3 = G3 + DOmega3_Dparam'*DOmega3_Dparam;                
            end
            for l3=l2:length(w); %loop computes analytical derivative of spectra and stacks them in a big matrix
                if ~Ident_Test.cumulants.fourth
                    break
                else
                    if rem(l3,100) == 0                    
                        msg = sprintf('l1=%d,l2=%d,l3=%d of %d each, with l1<l2<l3',l1,l2,l3,N+1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    elseif l1==length(w) && l2==l1 && l3==l2
                        msg = sprintf('l1=%d,l2=%d,l3=%d of %d each, with l1<l2<l3\n',l1,l2,l3,N+1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                    z3=exp(-1i*w(l3)); % Use Fourier transform for lag operator
                    [Hz3,DHz3_Dparam,DHz3cT_Dparam] = TransferFunction(z2,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam);
                    Hz2kronHz3 = kron(Hz2,Hz3);
                    DHz2kronHz3_Dparam = DerivXkronY(Hz2,DHz2_Dparam,Hz3,DHz3_Dparam);
                    [Hz1z2z3,DHz1z2z3_Dparam,DHz1z2z3cT_Dparam] = TransferFunction(z1*z2*z3,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam);
                    Hz1z2z3kronHz1 = kron(Hz1z2z3,Hz1);
                    DHz1z2z3kronHz1T_Dparam = DerivXkronY(Hz1z2z3',DHz1z2z3cT_Dparam,Hz1,DHz1cT_Dparam);
                    for j=1:nparam
                        dHz2kronHz3 = reshape(DHz2kronHz3_Dparam(:,j),size(Hz2kronHz3));
                        dGAMMA4full = reshape(DGAMMA4full_Dparam(:,j),size(GAMMA4full));
                        dHz1z2z3kronHz1T = reshape(DHz1z2z3kronHz1T_Dparam(:,j),size(Hz1z2z3kronHz1'));
                        dOmega4 = dHz2kronHz3*GAMMA4full*Hz1z2z3kronHz1' + Hz2kronHz3*dGAMMA4full*Hz1z2z3kronHz1' + Hz2kronHz3*GAMMA4full*dHz1z2z3kronHz1T;
                        DOmega4_Dparam(:,j) = dOmega4(:);
                    end
                    % Compute derivative of trispectrum
                    DOmega4_Dparam = (1/(8*pi^3))*DOmega4_Dparam;
                    %G4 = G4 + transpose(commutation(nd^2,nd^2)*DOmega4_Dparam)*DOmega4_Dparam;                    
                    G4 = G4 + DOmega4_Dparam'*DOmega4_Dparam;                    
                end                
            end            
            if rem(l2,100) == 0
                msg = sprintf('l1=%d,l2=%d of %d each, with l1<l2',l1,l2,N+1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
            elseif l1==length(w)&&l2==l1
                msg = sprintf('l1=%d,l2=%d of %d each, with l1<l2\n',l1,l2,N+1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
        if Ident_Test.cumulants.second && rem(l1,100) == 0
            msg = sprintf('l1=%d of %d',l1,N+1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
        elseif Ident_Test.cumulants.second && l1==length(w)
            msg = sprintf('l1=%d of %d\n',l1,N+1); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end;    
    
    G2 = 2*pi*G2./length(w); % Normalize G2 Matrix
    G3 = 4*pi^2*G3./(length(w)^2); % Normalize G3 Matrix
    G4 = 8*pi^3*G4./(length(w)^3); % Normalize G4 Matrix
else
    G2 = []; G3 = []; G4 = []; % If not Qu Tkachenko's test
end
Deriv.QT_G2 = G2;
Deriv.QT_G3 = G3;
Deriv.QT_G4 = G4;
Deriv.QT_G = G2+G3+G4;



function [H,DH_Dparam,DHcT_Dparam] = TransferFunction(z,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam)
    % Compute H and DH_Dparam and its conjugate(!) transpose.         
    zIminusA =  (z*speye(size(prun_A,1)) - prun_A);
    zIminusAinv = zIminusA\speye(size(prun_A,1));
    DzIminusA_Dparam = -Dprun_A_Dparam;
    DzIminusAinv_Dparam = kron(-(transpose(zIminusA)\speye(size(prun_A,1))),zIminusAinv)*DzIminusA_Dparam;
    H = SelectMat*(prun_D + prun_C*zIminusAinv*prun_B); % Transfer function
    DH_Dparam = kron(speye(size(prun_D,2)),SelectMat)*(Dprun_D_Dparam + DerivABCD(prun_C,Dprun_C_Dparam,zIminusAinv,DzIminusAinv_Dparam,prun_B,Dprun_B_Dparam));
    DHcT_Dparam = commutation(size(H))*conj(DH_Dparam); % conjugate transpose!
end%transferfunction end

end % main function end
