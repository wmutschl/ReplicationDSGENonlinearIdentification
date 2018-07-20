function [G,Omega] = gmatrix(Solut,Deriv,DSGE_Model,Ident_Test,Settings)

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


%% Construct Dprun_A_Dparam, Dprun_B_Dparam, Dprun_C_Dparam, Dprun_D_Dparam, Dprun_c_Dparam, Dprun_d_Dparam
hx = hv(Solut.ind.hx); hu = hv(Solut.ind.hu); hss=hSS(1:nx);
Hxx= hvv(Solut.ind.Hxx); Hxu=hvv(Solut.ind.Hxu); Hux=hvv(Solut.ind.Hux); Huu=hvv(Solut.ind.Huu);
Dhx_Dparam = Dhv_Dparam(Solut.ind.hx,:); Dhu_Dparam = Dhv_Dparam(Solut.ind.hu,:);
DHxx_Dparam=Dhvv_Dparam(Solut.ind.Hxx,:); 
DHxu_Dparam=Dhvv_Dparam(Solut.ind.Hxu,:); 
DHux_Dparam=Dhvv_Dparam(Solut.ind.Hux,:); 
DHuu_Dparam=Dhvv_Dparam(Solut.ind.Huu,:);
Dhss_Dparam=DhSS_Dparam(1:nx,:);
Dhu_kron_hu = DerivXkronY(hu,Dhu_Dparam,hu,Dhu_Dparam);
Dhu_kron_hx = DerivXkronY(hu,Dhu_Dparam,hx,Dhx_Dparam);
Dhx_kron_hu = DerivXkronY(hx,Dhx_Dparam,hu,Dhu_Dparam);    
gx = gv(Solut.ind.gx); gu = gv(Solut.ind.gu);
Gxx=gvv(Solut.ind.Gxx); Gxu=gvv(Solut.ind.Gxu); Gux=gvv(Solut.ind.Gux); Guu=gvv(Solut.ind.Guu);
Dgx_Dparam = Dgv_Dparam(Solut.ind.gx,:); Dgu_Dparam = Dgv_Dparam(Solut.ind.gu,:);
DGxx_Dparam=Dgvv_Dparam(Solut.ind.Gxx,:); 
DGxu_Dparam=Dgvv_Dparam(Solut.ind.Gxu,:); 
DGux_Dparam=Dgvv_Dparam(Solut.ind.Gux,:); 
DGuu_Dparam=Dgvv_Dparam(Solut.ind.Guu,:);

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
Deriv.Dprun_A_Dparam = Dprun_A_Dparam; Deriv.Dprun_B_Dparam = Dprun_B_Dparam; Deriv.Dprun_C_Dparam = Dprun_C_Dparam; Deriv.Dprun_D_Dparam = Dprun_D_Dparam;
Deriv.Dprun_c_Dparam = Dprun_c_Dparam; Deriv.Dprun_d_Dparam = Dprun_d_Dparam;
%% Analytical derivative of Expectation of observables (Ed)
[Ed,DEd_Dparam,E_xf_xf,DE_xf_xf_Dparam] = DerivExpectation(Solut,Deriv,DSGE_Model.numbers,'Analytical',Settings.approx);
Solut.Ed = Ed; Deriv.DEd_Dparam = DEd_Dparam;
Solut.E_xf_xf = E_xf_xf; Deriv.DE_xf_xf_Dparam = DE_xf_xf_Dparam;

%% Analytical Derivative of cumulants of Innovations xi_t=[u;kron(u,u)-vec(SIGU);kron(u,xf);kron(xf,u)]
[M2min,M3min,M4min,DM2min_Dparam,DM3min_Dparam,DM4min_Dparam] = ProdMom_inov(nu,nx,Sigma,DSigma_Dparam,E_xf_xf,DE_xf_xf_Dparam,dfstudt,Ddfstudt_Dparam,DSGE_Model,Ident_Test,Settings);

%% Save or load duplication matrix from file depending on speed setting
filename = ['./models/', DSGE_Model.shortname,'/',DSGE_Model.shortname,'_spec',num2str(DSGE_Model.spec),'_approx',];
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



%% Qu and Tkachenko's criteria
    % Construct G analytically
    fprintf('Compute Polyspectra for frequencies\n')
    % Create vector of Fourier frequencies for approximation of the integral    
    N=str2double(Ident_Test.options{2}); % Subintervalls    
    % Create Fourier frequencies
    w=2*pi*(-(N/2):1:(N/2))'/N;
    % Calculate zero-lag cumulants and its derivative
        % Auxiliary matrix that selects unique elements in xi_t
        Fxi = [speye(nu) spalloc(nu,nu*(nu+1)/2 + nu*nx,0);...
                spalloc(nu^2,nu,0) sparse(duplication(nu)) spalloc(nu^2,nu*nx,0);...
                spalloc(nu*nx,nu+nu*(nu+1)/2,0) speye(nu*nx); 
                spalloc(nx*nu,nu+nu*(nu+1)/2,0) sparse(commutation(nx,nu))]; 
       
        GAMMA2 = reshape(DPxi*GAMMA2min,size(Fxi,2),size(Fxi,2)); %Note GAMMA2 is a matrix, not a vector
        
        for i=1:nparam  
            dGAMMA2 = reshape(DPxi*DGAMMA2min_Dparam(:,i),size(Fxi,2),size(Fxi,2));
            DGAMMA2_Dparam(:,i) = dGAMMA2(:);
            clear dGAMMA2
        end
    GAMMA2full = Fxi*GAMMA2*transpose(Fxi);
    for j=1:nparam
        DGAMMA2full_Dparam(:,j) = vec(Fxi*reshape(DGAMMA2_Dparam(:,j),size(Fxi,2),size(Fxi,2))*transpose(Fxi));
    end
    G = zeros(nparam,nparam,length(w)); % Initialize objective function for power spectrum
    Omega = zeros(length(w),size(DSGE_Model.symbolic.SelectMat,1)^2); % Initialize objective function for power spectrum
    tic
    parfor l1=1:length(w); %loop computes analytical derivative of spectra
        z1=exp(-1i*w(l1)); % Use Fourier transform for lag operator
        [Hz1,DHz1_Dparam,DHz1cT_Dparam] = TransferFunction(z1,SelectMat,prun_A,prun_B,prun_C,prun_D,Dprun_A_Dparam,Dprun_B_Dparam,Dprun_C_Dparam,Dprun_D_Dparam);
        % Compute power spectrum
        Omega(l1,:) = (1/(2*pi))*transpose(vec(Hz1*GAMMA2full*Hz1'));
        % Compute derivative of power spectrum
        DOmega2_Dparam = (1/(2*pi))*DerivABCD(Hz1,DHz1_Dparam,GAMMA2full,DGAMMA2full_Dparam,Hz1',DHz1cT_Dparam);
        G(:,:,l1) = DOmega2_Dparam'*DOmega2_Dparam;
    end;
    toc
    G = 2*pi*sum(G,3)./length(w); % Normalize G2 Matrix    

