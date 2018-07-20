% Quadruplication Matrix as defined by
% Meijer (2005) - Matrix algebra for higher order moments. Linear Algebra and its Applications, 410,pp. 112–134
%
% Inputs:
%       p: size of vector
% Outputs:
%       QP: quadruplication matrix
%       QPinv: Moore-Penrose inverse of QP
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [QP,QPinv] = quadruplication(p)

QP = spalloc(p^4,p*(p+1)*(p+2)*(p+3)/24,p^4);
QPinv = spalloc(p*(p+1)*(p+2)*(p+3)/24,p*(p+1)*(p+2)*(p+3)/24,p^4);

for l=1:p
    for k=l:p
        for j=k:p
            for i=j:p                
                idx = perms([i j k l]);
                for r = 1:size(idx,1)
                    ii = idx(r,1); jj= idx(r,2); kk=idx(r,3); ll=idx(r,4);
                    n = ii + (jj-1)*p + (kk-1)*p^2 + (ll-1)*p^3;                    
                    m = mue(p,i,j,k,l);
                    QP(n,m)=1;
                    if i==j && j==k && k==l
                        QPinv(m,n)=1;
                    elseif i==j && j==k && k>l
                        QPinv(m,n)=1/4;
                    elseif i>j && j==k && k==l
                        QPinv(m,n)=1/4;
                    elseif i==j && j>k && k==l
                        QPinv(m,n) = 1/6;
                    elseif i>j && j>k && k==l
                        QPinv(m,n) = 1/12;
                    elseif i>j && j==k && k>l
                        QPinv(m,n) = 1/12;
                    elseif i==j && j>k && k>l
                        QPinv(m,n) = 1/12;
                    elseif i>j && j>k && k>l
                        QPinv(m,n) = 1/24;                    
                    end
               end
            end
        end
    end
end
%QPinv = (transpose(QP)*QP)\transpose(QP);

function m = mue(p,i,j,k,l)
     m = i + (j-1)*p + 1/2*(k-1)*p^2 + 1/6*(l-1)*p^3 - 1/2*j*(j-1) + 1/6*k*(k-1)*(k-2) - 1/24*l*(l-1)*(l-2)*(l-3) - 1/2*(k-1)^2*p + 1/6*(l-1)^3*p - 1/4*(l-1)*(l-2)*p^2 - 1/4*l*(l-1)*p + 1/6*(l-1)*p;
     m = round(m);
end


end