% Compare reconstructed indicator function for gmm 2d

rng('default');

N_list = [500; 1000; 1500; 2000; 2500; 3000; 3500;];
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
% closed form pdf
%         p = @(x,y) alpha(1)*mvnpdf([x y],mu_1,sigma_1) + ...
%                    alpha(2)*mvnpdf([x y],mu_2,sigma_2);
p = @(x,y) weights(1)*normpdf(x,mu_1(1),sqrt(sigma_1(1,1))).*normpdf(y,mu_1(2),sqrt(sigma_1(2,2))) + ...
           weights(2)*normpdf(x,mu_2(1),sqrt(sigma_2(1,1))).*normpdf(y,mu_2(2),sqrt(sigma_2(2,2)));

thresh = 10^-4;

for index1 = 1:reps
    tic 
    
    for index2 = 1:length(N_list)

        N = N_list(index2);
        
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


        % x = 0
        f = (X(:,1) > 0);
        f_gft = E_vec'*f;
        e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
        bw_index = find(e_dist > (1-thresh),1);
        bw_list(index1,index2,1) = lambda(bw_index);

        % x = -1
        f = (X(:,1) > -1);
        f_gft = E_vec'*f;
        e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
        bw_index = find(e_dist > (1-thresh),1);
        bw_list(index1,index2,2) = lambda(bw_index);

        % x = y^2 - 1
        f = (X(:,1) > X(:,2).^2 - 1);
        f_gft = E_vec'*f;
        e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
        bw_index = find(e_dist > (1-thresh),1);
        bw_list(index1,index2,3) = lambda(bw_index);

        % y = 0
        f = (X(:,2) > 0);
        f_gft = E_vec'*f;
        e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
        bw_index = find(e_dist > (1-thresh),1);
        bw_list(index1,index2,4) = lambda(bw_index);

        % x^2 + y^2 = 1
        f = (X(:,1).^2 + X(:,2).^2 <= 1);
        f_gft = E_vec'*f;
        e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
        bw_index = find(e_dist > (1-thresh),1);
        bw_list(index1,index2,5) = lambda(bw_index);

    %     progress_percent = index/reps*100;
    %     if mod(progress_percent,10) == 0
    %         fprintf('progress = %d%%\n', progress_percent);
    %     end
    end
    
    fprintf('Iteration: %d\n', index1);
    toc
end

mean_bw_list = squeeze(mean(bw_list));
std_bw_list = squeeze(std(bw_list));

rng('shuffle');

save('results/sep_exp_2.mat','bw_list','mean_bw_list','std_bw_list','p');

%% Plotting

addpath(genpath('plotting'));

font_size = 16;

figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.3]);

plot(N_list, std_bw_list, '--o', 'LineWidth',2);
xlabel('n','FontSize',font_size);
ylabel('std. dev.','FontSize',font_size);
set(gca,'FontSize',font_size);
xlim([0 4000]);
ylim([0 0.1]);

legend('\partial S_1','\partial S_2','\partial S_3','\partial S_4','\partial S_5');

export_fig('plots/sep_exp_2.pdf','-transparent');