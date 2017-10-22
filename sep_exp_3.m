% Compare reconstructed indicator function for gmm 2d

rng('default');

N = 2500;
sigma = 100/1000;

% Number of times experiment will be repeated
reps = 100;

tic        

% gaussians
mu_1 = [-1 0];
mu_2 = [1 0];
sigma_1 = 0.5^2*eye(2);
sigma_2 = 0.4^2*eye(2);
% mixing
weights = [0.4 0.6];
% closed form pdf
%         p = @(x,y) alpha(1)*mvnpdf([x y],mu_1,sigma_1) + ...
%                    alpha(2)*mvnpdf([x y],mu_2,sigma_2);
p = @(x,y) weights(1)*normpdf(x,mu_1(1),sqrt(sigma_1(1,1))).*normpdf(y,mu_1(2),sqrt(sigma_1(2,2))) + ...
           weights(2)*normpdf(x,mu_2(1),sqrt(sigma_2(1,1))).*normpdf(y,mu_2(2),sqrt(sigma_2(2,2)));

step = 0.05;
x = -3:step:3;
y = x;
[X,Y] = meshgrid(x,y);
Z = p(X(:),Y(:));
Z = reshape(Z,size(X));

num_cuts = 5;
sup_cut = [p(0,0) p(-1,0) p(-1,0) p(1,0) p(1,0)];

for j = 1:num_cuts
    lc_cut(j) = sum(Z(Z<=sup_cut(j))) *(step)^2;
end

label_percent_list = 0.05:0.05:0.95;

% Obtain best labeled sets using column-wise Gaussian elimination with
load('results/sep_exp_3_label_order.mat');

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

    for j = 1:length(label_percent_list)
%         K = round((round(sup_cut(1)*100) + label_percent_list(j))/100*N);
        K = round(label_percent_list(j)*N);
        L_set = label_order(1:K,index);
        U_set = setdiff(1:N,L_set);
        M = E_vec(:,1:K)*E_vec(L_set,1:K)^(-1);
        
        % x = 0
        f = (X(:,1) > 0);
        f_L = f(L_set);
        f_recon = (M*f_L > 0.5);
        error_list_1(j,index) = sum(abs(f(U_set)-f_recon(U_set)))/length(U_set);

        
        % x = -1
        f = (X(:,1) > -1);
        f_L = f(L_set);
        f_recon = (M*f_L > 0.5);
        error_list_2(j,index) = sum(abs(f(U_set)-f_recon(U_set)))/length(U_set);

        % x = y^2 - 1
        f = (X(:,1) > X(:,2).^2 - 1);
        f_L = f(L_set);
        f_recon = (M*f_L > 0.5);
        error_list_3(j,index) = sum(abs(f(U_set)-f_recon(U_set)))/length(U_set);
        
        % y = 0
        f = (X(:,2) > 0);
        f_L = f(L_set);
        f_recon = (M*f_L > 0.5);
        error_list_4(j,index) = sum(abs(f(U_set)-f_recon(U_set)))/length(U_set);

        % x^2 + y^2 = 1
        f = (X(:,1).^2 + X(:,2).^2 <= 1);
        f_L = f(L_set);
        f_recon = (M*f_L > 0.5);
        error_list_5(j,index) = sum(abs(f(U_set)-f_recon(U_set)))/length(U_set);
    end
    
    fprintf('Iteration: %d\n', index);

%     progress_percent = index/reps*100;
%     if mod(progress_percent,10) == 0
%         fprintf('progress = %d%%\n', progress_percent);
%     end

    toc
end

mean_error_list_1 = mean(error_list_1,2);
std_error_list_1 = std(error_list_1,[],2);

mean_error_list_2 = mean(error_list_2,2);
std_error_list_2 = std(error_list_2,[],2);

mean_error_list_3 = mean(error_list_3,2);
std_error_list_3 = std(error_list_3,[],2);

mean_error_list_4 = mean(error_list_4,2);
std_error_list_4 = std(error_list_4,[],2);

mean_error_list_5 = mean(error_list_5,2);
std_error_list_5 = std(error_list_5,[],2);

rng('shuffle');

save('results/sep_exp_3.mat','label_percent_list','lc_cut','p',...
'error_list_1','mean_error_list_1','std_error_list_1',...
'error_list_2','mean_error_list_2','std_error_list_2',...
'error_list_3','mean_error_list_3','std_error_list_3',...
'error_list_4','mean_error_list_4','std_error_list_4',...
'error_list_5','mean_error_list_5','std_error_list_5' );

%% Plotting

addpath('plotting');

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


plot(label_percent_list, mean_error_list_2, 'k-^', [lc_cut(2) lc_cut(2)], [0 0.1], 'r--','LineWidth',2);
ylabel('Mean error','FontSize',font_size);
xlabel('Fraction of labeled examples','FontSize',font_size);
set(gca,'FontSize',font_size);
ylim([0 0.1]);
set(gca,'YTick',0:0.02:0.1);
% export_fig('plots/sep_exp_3.pdf','-transparent');
