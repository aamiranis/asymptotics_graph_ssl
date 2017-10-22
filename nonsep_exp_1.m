% Compare reconstructed indicator function for gmm 2d

rng('default');

N = 2500;
sigma = 100/1000;

% Number of times experiment will be repeated
reps = 100;

tic        

x0 = 0.75;

% mixing
alpha = 3;
weights = [0.5 0.5];

p1 = @(x,y) (alpha/pi) * ( 1 - (x + x0).^2 - y.^2 ).^(alpha-1) .* ((x + x0).^2 + y.^2 <= 1);
p2 = @(x,y) (alpha/pi) * ( 1 - (x - x0).^2 - y.^2 ).^(alpha-1) .* ((x - x0).^2 + y.^2 <= 1);
p = @(x,y) weights(1)*p1(x,y) + weights(2)*p2(x,y);

% p = @(x,y) (alpha/pi) * ( 1 - x.^2 - y.^2 ).^(alpha-1) .* (x.^2 + y.^2 <= 1);

step = 0.025;
x = -2:step:2;
y = x;
[X_grid,Y_grid] = meshgrid(x,y);
Z = p(X_grid(:),Y_grid(:));
Z = reshape(Z,size(X_grid));
boundary = p1(X_grid,Y_grid) .* p2(X_grid,Y_grid) > 0;

sup_cut = max(max(Z .*  boundary));

bw_list = zeros(reps,1);

M = 20;
proxy_list = zeros(reps,M);
integral_values = zeros(M,1);

for m = 1:M
    integral_values(m) = ( weights(1) * weights(2) * ...
        sum( p1(X_grid(:),Y_grid(:)) .* p2(X_grid(:),Y_grid(:)) .* p(X_grid(:),Y_grid(:)).^(m-1) ) * (step)^2 )^(1/m);
end

thresh = 10^-4;

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
   
    f = logical(mix(:,1));
    
    f_gft = E_vec'*f;
    e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
    bw_index = find(e_dist > (1-thresh),1);
    bw_list(index) = lambda(bw_index);
 
    f_temp = f;
    for m = 1:M
        f_temp = L*f_temp;
        proxy_list(index,m) = (1/N * f'*f_temp)^(1/m);
    end
    
    
    fprintf('Iteration: %d\n', index);

%     progress_percent = index/reps*100;
%     if mod(progress_percent,10) == 0
%         fprintf('progress = %d%%\n', progress_percent);
%     end

%     save('results/exp_2.mat','p','sup_cut','bw_list','proxy_list','integral_values');
    
    toc
    
end

mean_bw_list = mean(bw_list);
std_bw_list = std(bw_list);

save('results/nonsep_exp_1.mat','p','sup_cut','bw_list','mean_bw_list','std_bw_list','proxy_list','integral_values');

%% Plotting

% addpath(genpath('~/matlab/plotting'));

font_size = 16;

figure1 = figure;

scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.45]);

step = 0.1;
x = -2:step:2;
y = -2:step:2;
[X_grid,Y_grid] = meshgrid(x,y);
Z = p(X_grid(:),Y_grid(:));
Z = reshape(Z,size(X_grid));
surf(X_grid,Y_grid,Z);
% shading interp;
view(345,10);

% zlim([0 0.6]);
set(gca,'ZTick',0:0.2:0.6);
set(gca,'YTick',-2:2:2);
set(gca,'XTick',-3:1:3);
set(gca, 'FontSize', font_size);
xlabel('x','FontSize', font_size);
ylabel('y','FontSize', font_size);
zlabel('p(x,y)','FontSize', font_size);

export_fig('plots/nonsep_density.pdf','-transparent');

close(figure1);