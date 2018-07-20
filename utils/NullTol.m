% Computes rank and analyzes null space of a matrix using the singular
% value decomposition for a specified tolerance level
%
% Inputs: X: matrix to be analyzed; tol: tolerance level
% Outputs: 
%       - rank_X: number of non-zero singular values which is the same as the
%       number of non-zero diagonal elements in S
%       - Null_X: Null space of X, since the right-singular vectors (columns of V) corresponding to vanishing singular values of X span the null space of X
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [rank_X,Null_X,tol] = NullTol(X,tol)

    if ~isnumeric(tol)
        tol = max(size(X))*eps(norm(X));
    else
        tol = tol;
    end
    [U,S,V] = svd(X,0);
    if size(S,2) == 1
        rank_X = sum(S > tol);
        Null_X = 0;    
    else
        rank_X = sum(diag(S) > tol);
        Null_X = V(:,find(abs(diag(S))<tol));    
    end    
end %NullTol end