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

step = 0.025;
x = -2:step:2;
y = x;
[X_grid,Y_grid] = meshgrid(x,y);
Z = p(X_grid(:),Y_grid(:));
Z = reshape(Z,size(X_grid));
boundary = p1(X_grid,Y_grid) .* p2(X_grid,Y_grid) > 0;

sup_cut = max(max(Z .*  boundary));

lc_cut = sum(Z(Z<=sup_cut)) *(step)^2;

label_percent_list = 0.05:0.05:0.95;

% Obtain best labeled sets using column-wise Gaussian elimination with
load('results/nonsep_exp_3_label_order.mat');

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

    for j = 1:length(label_percent_list)
%         K = round((round(sup_cut(1)*100) + label_percent_list(j))/100*N);
        K = round(label_percent_list(j)*N);
        L_set = label_order(1:K,index);
        U_set = setdiff(1:N,L_set);
        M = E_vec(:,1:K)*E_vec(L_set,1:K)^(-1);
        
        f = logical(mix(:,1));
        f_L = f(L_set);
        f_recon = (M*f_L > 0.5);
        error_list(j,index) = sum(abs(f(U_set)-f_recon(U_set)))/length(U_set);

%     progress_percent = index/reps*100;
%     if mod(progress_percent,10) == 0
%         fprintf('progress = %d%%\n', progress_percent);
%     end
    end

    fprintf('Iteration: %d\n', index);
    
    toc
end

mean_error_list = mean(error_list,2);
std_error_list = std(error_list,[],2);

rng('shuffle');

save('results/nonsep_exp_3.mat','label_percent_list','lc_cut','sup_cut','p',...
'error_list','mean_error_list','std_error_list' );


%% Plotting

addpath('~/matlab/plotting');

font_size = 16;

figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.3]);

% plot(label_percent_list, mean_error_list_1, 'b-o', [lc_cut(1) lc_cut(1)], [0 0.01], 'r--','LineWidth',2);
% ylabel('Mean error','FontSize',font_size);
% xlabel('Fraction of labeled examples','FontSize',font_size);
% ylim([0 0.01]);
% set(gca,'YTick',0:0.002:0.01);
% set(gca,'FontSize',font_size);
% export_fig('plots_aistats/exp3_1.pdf','-transparent');

% plot(label_percent_list, mean_error_list_2, 'k-^',label_percent_list, mean_error_list_3, 'b-o', [lc_cut(2) lc_cut(2)], [0 0.1], 'r--','LineWidth',2);
% ylabel('Mean error','FontSize',font_size);
% xlabel('Fraction of labeled examples','FontSize',font_size);
% set(gca,'FontSize',font_size);
% ylim([0 0.1]);
% set(gca,'YTick',0:0.02:0.1);
% l1 = legend('\partialS_2','\partialS_3');
% set(l1,'Position',[0.70 0.71 0.20 0.07]);
% export_fig('plots_aistats/exp3_2.pdf','-transparent');

% plot(label_percent_list, mean_error_list_4, 'k-^',label_percent_list, mean_error_list_5, 'b-o', [lc_cut(4) lc_cut(4)], [0 0.4], 'r--','LineWidth',2);
% ylabel('Mean error','FontSize',font_size);
% xlabel('Fraction of labeled examples','FontSize',font_size);
% set(gca,'FontSize',font_size);
% ylim([0 0.4]);
% xlim([0 1.01]);
% set(gca,'YTick',0:0.1:0.4);
% l1 = legend('\partialS_4','\partialS_5');
% set(l1,'Position',[0.60 0.71 0.20 0.07]);
% export_fig('plots_aistats/exp3_3.pdf','-transparent');


plot(label_percent_list, mean_error_list, 'k-^', [lc_cut lc_cut], [0 0.1], 'r--','LineWidth',2);
ylabel('Mean error','FontSize',font_size);
xlabel('Fraction of labeled examples','FontSize',font_size);
set(gca,'FontSize',font_size);
ylim([0 0.1]);
set(gca,'YTick',0:0.02:0.1);
% export_fig('plots/nonsep_exp_3.pdf','-transparent');
