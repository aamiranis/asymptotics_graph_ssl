% compute pairwise distances
perm1 = zeros(N*(N-1)/2,1);
perm2 = zeros(N*(N-1)/2,1);
perm3 = zeros(N*(N-1)/2,1);
k = 1;
for i = 1:N-1
    for j = i+1:N
        perm1(k) = i;
        perm2(k) = j;
        perm3(k) = N*(i-1) + j;
        k = k+1;
    end
end

d = size(X,2);
if (d > 1)
    temp = sum((X(perm1,:) - X(perm2,:)).^2,2);
else
    temp = (X(perm1) - X(perm2)).^2;
end
temp = 1/(2*pi*sigma^2)^(d/2) * exp(-temp/(2*sigma^2));

% Complete upper triangular part
% do this after adjacency is computed
L = zeros(N);
L(perm3) = temp;
% clear temp;

% creating adjacency
L = L + L';

% flag_sparse = false;
% if (nnz(L) < N*N/2)
%     L = sparse(L);
%     flag_sparse = true;
% end

temp = sum(L,2);

% % creating laplacian
% if (flag_sparse)
%     L = 1/N*(spdiags(temp,0,N,N) - L);
% else
%     L = 1/N*(diag(temp) - L);
% end

L = 1/N*(diag(temp) - L);

% % Full matrices
% temp = sum(L,2);
% L = 1/N*(diag(temp) - L);

% % If eigenvalue decomposition is required
% Lam = sort(eig(L));