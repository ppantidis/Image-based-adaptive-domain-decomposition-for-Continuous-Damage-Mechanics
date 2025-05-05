function[Khh_inv,CA_inverse,CA_inverseB,lenA,lenD] = func_getSchurBlocks(K_F, Domain_Vector,listofnodes_ebc)

count = 0;        % for checking what to exclude from K matrix
len_noDamage = 0; % size of K matrix where no damage is going to occur
j = 1;

while ~ Domain_Vector(j).DomainDamage
    len_noDamage = len_noDamage + length(Domain_Vector(j).K);
    % this is done under the assumption that the domain that doesn't have damage comes before the one with damage
    j = j + 1;
end

for k = 1:size(listofnodes_ebc,1)
    if listofnodes_ebc(k) <= len_noDamage
        count = count + 1;  % cos K_F does not include nodes on essential boundary
                            % so keep track of them and then take them out
    end
end

% Get length of blocks and use the information to isolate the length of the blocks
lenA = length(Domain_Vector(1).K) - count;
lenD = length(K_F) - lenA;

% Compute the inverse of the healthy part of K
Khh_inv = (K_F(1:lenA, 1:lenA)) \ eye(size(K_F(1:lenA, 1:lenA)));
% Khh_inv = sparse(K_F(1:lenA, 1:lenA)) \ eye(size(K_F(1:lenA, 1:lenA)));

% Compute additional useful (and constant) matrices
CA_inverse  = K_F(lenA+1:lenA+lenD, 1:lenA) * Khh_inv; % 80 x 2812
CA_inverseB = CA_inverse * K_F(1:lenA, lenA+1:lenA+lenD);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using LU decomposition
% [L,U,P] = lu(K_F(1:lenA,1:lenA));

% % Using mldivide alteration
% CA_inverse =  (K_F(1:lenA,1:lenA)' \ K_F(lenA+1:lenA+lenD, 1:lenA)')';

% Using LU decomposition and forward/backward substitution
% Khh_inv = zeros(size(K_F(1:lenA,1:lenA)));
% I = eye(size(K_F(1:lenA,1:lenA)));
% for j = 1:length(I)
%     Khh_inv(:,j) = lup_solve(L,U,P,I(:,j));
% end
% CA_inverse = K_F(lenA+1:lenA+lenD, 1:lenA) * Khh_inv;

% Using mldivide over identity matrix
% Khh_inv = K_F(1:lenA, 1:lenA)\eye(size(K_F(1:lenA, 1:lenA)));
% CA_inverse = K_F(lenA+1:lenA+lenD, 1:lenA) * Khh_inv;

% % Using direct inversion
% CA_inverse  = K_F(lenA+1:lenA+lenD, 1:lenA) * inv(K_F(1:lenA, 1:lenA));

% tic
% Khh_inv     = K_F(1:lenA, 1:lenA) \ eye(size(K_F(1:lenA, 1:lenA)));
% CA_inverse  = K_F(lenA+1:lenA+lenD, 1:lenA) * Khh_inv;
% toc




%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     % ---------------------------------------------------------------------
%     function [x] = back_sub(U,b)
%         n = length(b);
%         x(n,1) = b(n) / U(n,n);
%         for i = n-1:-1:1
%             x(i,1) = (b(i) - U(i,i+1:n) * x(i+1:n,1)) ./ U(i,i);
%         end
%     end
% 
%     % ---------------------------------------------------------------------
%     function [x] = forward_sub(L,b)
%         n = length(L);
%         x = zeros(n,1);
%         for i = 1:n
%             x(i) = (b(i) - L(i,1:i-1)*x(1:i-1) / L(i,i));
%         end
%     end
% 
%     % ---------------------------------------------------------------------
%     function [x] = lu_solve(L,U,b)
%         y = forward_sub(L,b);
%         x = back_sub(U,y);
%     end
%     
%     % ---------------------------------------------------------------------
%     function [x] = lup_solve(L,U,P,b)
%         z = P * b;
%         x = lu_solve(L,U,z);
%     end

end






