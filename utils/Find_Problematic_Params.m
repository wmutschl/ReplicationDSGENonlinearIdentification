% Finds problematic sets of paramters given an objective function either
% via a brute-force approach (consider all possible subsets) or via the
% nullspace
%
% Input:    A:          Objective Matrix
%           Null_A:     Null space of objective matrix
%           tol:        tolerance level for rank calculations
%           parnames:   symbolic names of parameters
%           nparam:     number of parameters
%           proc:       procedure to calculate problematic parameters (bruteforce vs. Nullspace)
%           test:       which test to consider
%
% Output:   problpars: structure containing sets of problematic parameter
%
% Calls:    NullTol
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function problpars = Find_Problematic_Params(A,Null_A,tol,parnames,nparam,proc,test)
switch test
    case 'Iskrev'
        n_fixed_rank =0; %no additional rank number
        j=1; ElimPar=[]; %initialize, start with one-element sets
    case 'QuTkachenko'        
        n_fixed_rank =0;
        j=2; %Start with 2-element subsets
        % Check 1-element subsets
        ElimPar = find(abs(diag(A))<tol);
        if isempty(ElimPar) == 0
            problpars{1} = parnames(ElimPar);
        end;        
end

switch proc
    %% Better and more accurate, but slower procedure, since it checks all possible subsets
    case '1' % Brute Force procedure
        ElimParj=[];
        while j<=nparam % Check j-element subsets
            ElimPar = sort([ElimPar;ElimParj]);
            indexj=nchoosek(1:nparam,j);  %all possible subsets of j elements
            % Remove already problematic parameters
            for c=1:j
                for n=1:length(ElimPar)
                    indexj(indexj(:,c)==ElimPar(n),:)=[];
                end
            end
            % Initialize blanks
            rankj=zeros(size(indexj,1),1);
            if strcmp(test,'QuTkachenko')
                A_j=zeros(j,j,size(indexj,1));
            end
            %Check rank criteria
            for k=1:size(indexj,1)                
                switch test
                    case 'Iskrev'
                        A_j = A(:,indexj(k,:));
                    case 'QuTkachenko'
                        for ii=1:j
                            for jj=1:j
                            %extract the j-element Gs from the G matrix
                            A_j(ii,jj,k) = A(indexj(k,ii),indexj(k,jj)); 
                            end
                        end
                end
                % Compute rank with imposed tol                
                if strcmp(test,'QuTkachenko')
                    [rankj(k,1),~,~] = NullTol(A_j(:,:,k),tol);
                else
                    [rankj(k,1),~,~] = NullTol(A_j,tol);
                end
            end
            ProblPar = indexj(rankj < j+n_fixed_rank,:); %Compare rank condition for all possible subsets
            if isempty(ProblPar) == 0
                problpars{j} = parnames(ProblPar);
                if size(ProblPar,1)==1 %if there is only one set, put it into a row
                    problpars{j} = transpose(problpars{j});
                end
            end
            % Add parameters to already problematic parameters    
            ElimParj = unique(ProblPar(:));
            j=j+1;
        end
    %% This procedure is faster, but not as accurate 
    case '2' %Nullspace procedure
        ind_par = abs(Null_A)>tol;
        ind_problpar = transpose(unique(ind_par','rows'));        
        for i=1:nparam
            indx_col = find(sum(ind_problpar)==i);
            ProblPar={};
            for col = indx_col
                ProblPar = [ProblPar;parnames(ind_problpar(:,col))'];
            end
            problpars{i} = ProblPar;
        end
end

end %function end