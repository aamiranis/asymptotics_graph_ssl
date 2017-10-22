% Compare reconstructed indicator function for gmm 2d

rng('default');

N_list = [500; 1000; 1500; 2000; 2500; 3000; 3500;];
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

thresh = 10^-4;

for index1 = 1:reps
    tic 
    
    for index2 = 1:length(N_list)

        N = N_list(index2);
        
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

        % x = 0
        f = logical(mix(:,1));
        f_gft = E_vec'*f;
        e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
        bw_index = find(e_dist > (1-thresh),1);
        bw_list(index1,index2) = lambda(bw_index);

    %     progress_percent = index/reps*100;
    %     if mod(progress_percent,10) == 0
    %         fprintf('progress = %d%%\n', progress_percent);
    %     end
    end
    
    fprintf('Iteration: %d\n', index1);
    toc
end

mean_bw_list = mean(bw_list);
std_bw_list = std(bw_list);

rng('shuffle');

save('results/nonsep_exp_2.mat','N_list','bw_list','mean_bw_list','std_bw_list','p');

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

export_fig('plots/nonsep_exp_2.pdf','-transparent');