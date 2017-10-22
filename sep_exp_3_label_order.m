% Compare reconstructed indicator function for gmm 2d

rng('default');

N = 2500;
sigma = 100/1000;

% Number of times experiment will be repeated
reps = 100;

% gaussians
mu_1 = [-1 0];
mu_2 = [1 0];
sigma_1 = 0.5^2*eye(2);
sigma_2 = 0.4^2*eye(2);
% mixing
weights = [0.4 0.6];

label_order = zeros(N,reps);

for index = 1:reps

    tic 
    
    mix = mnrnd(1,weights,N);
    X = zeros(N,2);
    % generate dimension 1
    X(logical(mix(:,1)),1) = normrnd(mu_1(1),sqrt(sigma_1(1,1)),sum(mix(:,1)),1);
    X(logical(mix(:,2)),1) = normrnd(mu_2(1),sqrt(sigma_2(1,1)),sum(mix(:,2)),1);
    % generate dimension 2
    X(logical(mix(:,1)),2) = normrnd(mu_1(2),sqrt(sigma_1(2,2)),sum(mix(:,1)),1);
    X(logical(mix(:,2)),2) = normrnd(mu_2(2),sqrt(sigma_2(2,2)),sum(mix(:,2)),1);

    compute_L;

    % Eigenvalue decomposition
    [E_vec,E_val] = eig(L);
    [lambda,perm] = sort(abs(diag(E_val)));
    E_val = E_val(perm,perm);    
    E_vec = E_vec(:,perm);
    
    % Obtain best labeled sets using column-wise Gaussian elimination with
    % row-pivoting on eigenvector matrix
    % - column-wise needed because we to find linearly indep. rows
    % - row-pivoting needed, because we don't want to consult future
    %   eigenvectors, as we generally need it for U_{:,K}.
    label_order(:,index) = col_ge_row_pivot(E_vec);
    
    fprintf('Iteration: %d\n', index);

%     progress_percent = index/reps*100;
%     if mod(progress_percent,10) == 0
%         fprintf('progress = %d%%\n', progress_percent);
%     end

    toc
end

rng('shuffle');

save('results/sep_exp_3_label_order.mat','label_order');