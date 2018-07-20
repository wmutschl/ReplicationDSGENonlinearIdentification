% Triplication Matrix as defined by
% Meijer (2005) - Matrix algebra for higher order moments. Linear Algebra and its Applications, 410,pp. 112–134
%
% Inputs:
%       p: size of vector
% Outputs:
%       TP: triplication matrix
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function TP = triplication(p)

TP = zeros(p^3,p*(p+1)*(p+2)/6);

for k=1:p
    for j=k:p
        for i=j:p            
            idx = unique(perms([i j k]),'rows');
            for r = 1:size(idx,1)
                ii = idx(r,1); jj= idx(r,2); kk=idx(r,3);
                n = ii + (jj-1)*p + (kk-1)*p^2;                                    
                m = mue(p,i,j,k);
                TP(n,m)=1;                
            end
            n=n+1;
        end
    end
end
%TP_MP = (transpose(TP)*TP)\transpose(TP);

function m = mue(p,i,j,k)
    m = i+(j-1)*p + 1/2*(k-1)*p^2 - 1/2*j*(j-1) + 1/6*k*(k-1)*(k-2) - 1/2*(k-1)^2*p;
    m = round(m);
end

end