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