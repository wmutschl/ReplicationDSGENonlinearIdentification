% This function
% (1) evaluates all model objects and derivatives at steady-state values numerically
% (2) solves the model using the algorithm described in Gomme & Klein (2011) - Second order approximation of dynamic models without the use of tensors, in: Journal of Economic Dynamics and Control 35.4, pp. 604-615,
% (3) Sets up the pruned state-space representation
% (4) writes all solution matrices as sparse for further evaluation.
% (5) uses index matrices to keep track of terms belonging to states, shocks and controls for solution matrices
%
% Inputs:
%       param_estim: vector with numerical values for all parameters
%       DSGE_model: structure containing DSGE model independent of names
%       approx: order of approximation
%       Deriv_type: Analytical or numerical derivatives
%
% Calls:
%       numeval: Numerically evaluate analytical model objects that were written into files
%       SolveModel: Solve DSGE model up to order 2
%
% Outputs:
%       Solut: Structure containing sparse numerical evaluated solution matrices, expectation, autocovariogram, and spectrum 
%       Deriv: Structure containing sparse numerical evaluated analytical derivative of solution matrices, expectation, autocovariogram, and spectrum w.r.t. param_identif (up to order 2)
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [Solut,Derivs] = EvaluateSparse(param_estim,DSGE_Model,approx,Deriv_type)
%% Evaluate all objects and derivatives at steady-state values numerically
[ngra,nf,nhes,nSigma,netatilde,nsig,nSelectMat,nSS,ndfstudt,nd,nx,nv,ny,nu,...
    nDgra_Dparam,nDhes_Dparam,nDSigma_Dparam,nDetatilde_Dparam,nDsig_Dparam,nDsigetatilde_Dparam,nDSS_Dparam,nDdfstudt_Dparam]...
    = numeval(param_estim,DSGE_Model,approx,Deriv_type);
%% Solve the model using the procedure by Gomme & Klein (2011)
[gv,gvv,gSS,hv,hvv,hSS,SolM,SolN,SolQ,SolR,SolS,SolU,~] = SolveModel(ngra,nhes,netatilde*transpose(netatilde),nv,ny,approx);

%% Separate matrices into states and shocks
% Use notation and formulas of Andreasen et al (2014)

ind.hv = reshape(1:nv^2,nv,nv);
ind.hvv = reshape(1:nv^3,[nv^2 nv]);
ind.hvv_tensor=permute(reshape(ind.hvv,[nv nv nv]),[2 1 3]);
ind.hx = ind.hv(1:nx,1:nx);                             
ind.hu = ind.hv(1:nx,(nx+1):end);    
ind.hxx = ind.hvv_tensor(1:nx,1:nx,1:nx);
ind.huu = ind.hvv_tensor(1:nx,(nx+1):end,(nx+1):end);
ind.hxu = ind.hvv_tensor(1:nx,(nx+1):end,1:nx);
ind.hux = ind.hvv_tensor(1:nx,1:nx,(nx+1):end);
ind.Hxx = reshape(ind.hxx,nx,nx*nx);                
ind.Hxu = reshape(ind.hxu,nx,nx*nu);                
ind.Hux = reshape(ind.hux,nx,nu*nx);                
ind.Huu = reshape(ind.huu,nx,nu*nu);

ind.gv = reshape(1:ny*nv,ny,nv);
ind.gvv = reshape(1:ny*nv^2,[ny*nv nv]);
ind.gvv_tensor=permute(reshape(ind.gvv,[nv ny nv]),[2 1 3]);
ind.gx = ind.gv(:,1:nx);
ind.gu = ind.gv(:,(nx+1):end);    
ind.gxx = ind.gvv_tensor(1:ny,1:nx,1:nx);
ind.guu = ind.gvv_tensor(1:ny,(nx+1):end,(nx+1):end);
ind.gux = ind.gvv_tensor(1:ny,1:nx,(nx+1):end);
ind.gxu = ind.gvv_tensor(1:ny,(nx+1):end,1:nx);  
ind.Gxx = reshape(ind.gxx,ny,nx*nx);
ind.Gxu = reshape(ind.gxu,ny,nx*nu);
ind.Gux = reshape(ind.gux,ny,nu*nx);
ind.Guu = reshape(ind.guu,ny,nu*nu);

hx = hv(ind.hx); hu = hv(ind.hu); hss = hSS(1:nx);
Hxx = hvv(ind.Hxx); Hxu = hvv(ind.Hxu); Hux = hvv(ind.Hux); Huu = hvv(ind.Huu);
gx = gv(ind.gx); gu = gv(ind.gu); gss = gSS;
Gxx = gvv(ind.Gxx); Gxu = gvv(ind.Gxu); Gux = gvv(ind.Gux); Guu = gvv(ind.Guu);
    
   
%% The pruned state-space representation
prun_A = [hx                ,zeros(nx,nx)       ,zeros(nx,nx^2);
          zeros(nx,nx)      ,hx                 ,0.5*Hxx;
          zeros(nx*nx,nx)   ,zeros(nx*nx,nx)    ,kron(hx,hx)  ];
prun_B = [hu zeros(nx,nu^2) zeros(nx,nu*nx) zeros(nx,nu*nx);...
          zeros(nx,nu) 0.5*Huu 0.5*Hux 0.5*Hxu;...
          zeros(nx^2,nu) kron(hu,hu) kron(hu,hx) kron(hx,hu)];
prun_C = [gx gx 0.5*Gxx];
prun_D = [gu 0.5*Guu 0.5*Gux 0.5*Gxu];
prun_c = [zeros(nx,1);
          0.5*hss*nsig^2 + 0.5*Huu*nSigma(:);
          kron(hu,hu)*nSigma(:)];
prun_d = 0.5*gss*nsig^2 + 0.5*Guu*nSigma(:);

%% Make system matrices sparse
Solut.f = sparse(nf); 
Solut.SelectMat=sparse(nSelectMat);
Solut.gra = sparse(ngra);                    Derivs.Dgra_Dparam = sparse(nDgra_Dparam); 
Solut.hes=sparse(nhes);                      Derivs.Dhes_Dparam=sparse(nDhes_Dparam); 
Solut.Sigma=sparse(nSigma);                  Derivs.DSigma_Dparam=sparse(nDSigma_Dparam);
Solut.etatilde=sparse(netatilde);            Derivs.Detatilde_Dparam=sparse(nDetatilde_Dparam);
Solut.SS = sparse(nSS);                      Derivs.DSS_Dparam = sparse(nDSS_Dparam);
Solut.sig = sparse(nsig);                    Derivs.Dsig_Dparam = sparse(nDsig_Dparam);
Solut.dfstudt = sparse(ndfstudt);            Derivs.Ddfstudt_Dparam = sparse(nDdfstudt_Dparam);
                                             Derivs.Dsigetatilde_Dparam=sparse(nDsigetatilde_Dparam);
Solut.gv = sparse(gv); Solut.gvv=sparse(gvv); Solut.gSS=sparse(gSS);
Solut.hv = sparse(hv); Solut.hvv=sparse(hvv); Solut.hSS=sparse(hSS);
Solut.ind = ind;
Solut.solM=sparse(SolM); Solut.solN=sparse(SolN); Solut.solQ=sparse(SolQ); Solut.solR=sparse(SolR);  Solut.solS=sparse(SolS);  Solut.solU=sparse(SolU);
Solut.prun_A = sparse(prun_A); Solut.prun_B = sparse(prun_B); Solut.prun_C = sparse(prun_C); Solut.prun_D = sparse(prun_D);
Solut.prun_c = sparse(prun_c);  Solut.prun_d = sparse(prun_d);
end