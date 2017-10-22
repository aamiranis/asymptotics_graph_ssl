% Compare reconstructed indicator function for gmm 2d

rng('default');

N = 2500;
sigma = 100/1000;

% Number of times experiment will be repeated
reps = 100;

x0 = 0.75;

% mixing
alpha = 3;
weights = [0.5 0.5];

p1 = @(x,y) (alpha/pi) * ( 1 - (x + x0).^2 - y.^2 ).^(alpha-1) .* ((x + x0).^2 + y.^2 <= 1);
p2 = @(x,y) (alpha/pi) * ( 1 - (x - x0).^2 - y.^2 ).^(alpha-1) .* ((x - x0).^2 + y.^2 <= 1);
p = @(x,y) weights(1)*p1(x,y) + weights(2)*p2(x,y);

label_order = zeros(N,reps);

for index = 1:reps

    tic
    
    mix = mnrnd(1,weights,N);
    X = zeros(N,2);
    
    N1 = sum(mix(:,1));
    R1 = sqrt( 1 - rand(N1,1).^(1/alpha) );
    T1 = rand(N1,1);
    
    N2 = sum(mix(:,2));
    R2 = sqrt( 1 - rand(N2,1).^(1/alpha) );
    T2 = rand(N2,1);
    
    % generate dimension 1
    X(logical(mix(:,1)),1) = R1 .* cos(2*pi*T1) - x0;
    X(logical(mix(:,2)),1) = R2 .* cos(2*pi*T2) + x0;
    
    % generate dimension 2
    X(logical(mix(:,1)),2) = R1 .* sin(2*pi*T1);
    X(logical(mix(:,2)),2) = R2 .* sin(2*pi*T2);

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

save('results/nonsep_exp_3_label_order.mat','label_order');