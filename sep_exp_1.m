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

% x_cut = [-1 0 1];
% num_cuts = length(x_cut);
% sup_cut = [p(-1,0) p(0,0) p(1,0)];

num_cuts = 5;
sup_cut = [p(0,0) p(-1,0) p(-1,0) p(1,0) p(1,0)];

% for j = 1:num_cuts
%     lc_cut(j) = sum(Z(Z<=sup_cut(j))) *(step)^2;
% end

bw_list = zeros(num_cuts,reps);

thresh = 10^-4;

for index = 1:reps

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
    
%     for j = 1:num_cuts
%         f = (X(:,1) > x_cut(j));            
%         f_gft = E_vec'*f;
%         e_dist = cumsum(f_gft.^2/sum(f_gft.^2));       
% %         lc = find(abs(f_gft) > 10^-2*max(abs(f_gft)),1,'last')/N;
%         lc = find(e_dist > (1-thresh),1)/N;       
%         lc_list(j,index) = lc;
%         f_gft_list(:,index,j) = f_gft;       
%     end

    % x = 0
    f = (X(:,1) > 0);
    f_gft = E_vec'*f;
    e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
    bw_index = find(e_dist > (1-thresh),1);
    bw_list(1,index) = lambda(bw_index);

    % x = -1
    f = (X(:,1) > -1);
    f_gft = E_vec'*f;
    e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
    bw_index = find(e_dist > (1-thresh),1);
    bw_list(2,index) = lambda(bw_index);
    
    % x = y^2 - 1
    f = (X(:,1) > X(:,2).^2 - 1);
    f_gft = E_vec'*f;
    e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
    bw_index = find(e_dist > (1-thresh),1);
    bw_list(3,index) = lambda(bw_index);
 
    % y = 0
    f = (X(:,2) > 0);
    f_gft = E_vec'*f;
    e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
    bw_index = find(e_dist > (1-thresh),1);
    bw_list(4,index) = lambda(bw_index);
    
    % x^2 + y^2 = 1
    f = (X(:,1).^2 + X(:,2).^2 <= 1);
    f_gft = E_vec'*f;
    e_dist = cumsum(f_gft.^2/sum(f_gft.^2));
    bw_index = find(e_dist > (1-thresh),1);
    bw_list(5,index) = lambda(bw_index);
    
    fprintf('Iteration: %d\n', index);

%     progress_percent = index/reps*100;
%     if mod(progress_percent,10) == 0
%         fprintf('progress = %d%%\n', progress_percent);
%     end
    
end

mean_bw_list = mean(bw_list,2);
std_bw_list = std(bw_list,[],2);

rng('shuffle');

save('results/sep_exp_1.mat','sup_cut','bw_list','mean_bw_list','std_bw_list','p');

%% Plotting

addpath(genpath('~/matlab/plotting'));

font_size = 16;

figure1 = figure;

scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.45]);

step = 0.1;
x = -3:step:3;
y = -2:step:2;
[X,Y] = meshgrid(x,y);
Z = p(X(:),Y(:));
Z = reshape(Z,size(X));
surf(X,Y,Z);
% shading interp;
view(345,10);

zlim([0 0.6]);
set(gca,'ZTick',0:0.2:0.6);
set(gca,'YTick',-2:2:2);
set(gca,'XTick',-3:1:3);
set(gca, 'FontSize', font_size);
xlabel('x','FontSize', font_size);
ylabel('y','FontSize', font_size);
zlabel('p(x,y)','FontSize', font_size);
export_fig('plots/sep_density.pdf','-transparent');

close(figure1);


figure2 = figure;
% scr = get(0,'ScreenSize');
% set(gcf,'PaperPositionMode','auto');
% set(figure2, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.45]);

contour(X,Y,Z,20,'LineWidth',1.5);
% colormap(hot);
xlim([-2 2]);
ylim([-2 2]);
axis equal;
set(gca,'FontSize',font_size);
set(gca,'XTick',-2:2);
set(gca,'YTick',-2:2);
xlabel('x', 'FontSize', font_size);
ylabel('y', 'FontSize', font_size);
% colorbar; caxis([0 0.6]);
hold all;

yy = [-2 2];
xx = [0 0];
plot(xx,yy,'k--','LineWidth',2);

yy = [-2 2];
xx = [-1 -1];
plot(xx,yy,'k--','LineWidth',2);

yy = [-2:0.1:2];
xx = yy.^2 - 1;
plot(xx,yy,'k--','LineWidth',2);

yy = [0 0];
xx = [-3 3];
plot(xx,yy,'k--','LineWidth',2);

yy = sin(0:(2*pi/24):2*pi);
xx = cos(0:(2*pi/24):2*pi);
plot(xx,yy,'k--','LineWidth',2);

% Create textbox
annotation(figure2,'textbox',[0.52 0.85 0.018 0.048],'String',{'\partialS_1'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure2,'textbox',[0.28 0.85 0.018 0.048],'String',{'\partialS_2'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure2,'textbox',[0.75 0.79 0.018 0.048],'String',{'\partialS_3'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure2,'textbox',[0.82 0.56 0.018 0.048],'String',{'\partialS_4'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure2,'textbox',[0.42 0.76 0.018 0.048],'String',{'\partialS_5'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');

export_fig('plots/sep_boundaries.pdf','-transparent');


figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.45]);

ax1 = bar([mean_bw_list+std_bw_list/2 sup_cut'],'LineWidth',1);
set(ax1(1),'FaceColor',[0.6 0.6 0.6]);
set(ax1(2),'FaceColor',[1 0.5 0.5]);
ylim([0 1.2]);
set(gca,'FontSize',font_size);
set(gca,'XTickLabel',{'','','','',''});
annotation(figure1,'textbox',[0.185 0.055 0.018 0.048],'String',{'\partialS_1'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.34 0.055 0.018 0.048],'String',{'\partialS_2'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.495 0.055 0.018 0.048],'String',{'\partialS_3'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.65 0.055 0.018 0.048],'String',{'\partialS_4'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.805 0.055 0.018 0.048],'String',{'\partialS_5'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
set(gca,'YTick',0:0.2:1);
hold all;

ax2 = bar([mean_bw_list-std_bw_list/2 sup_cut'],'LineWidth',1);
set(ax2(1),'FaceColor',[0.9 0.9 0.9]);
set(ax2(2),'FaceColor',[1 0.5 0.5]);

ax3 = bar([mean_bw_list sup_cut'],'LineWidth',1);
set(ax3(1),'FaceColor','none');
set(ax3(2),'FaceColor',[1 0.5 0.5]);

ylabel('Bandwidth');

legend1 = legend([ax2(1) ax2(2)],{'Empirical', 'Theoretical'});
set(legend1,...
    'Position',[0.15 0.77 0.28 0.08]);

export_fig('plots/sep_exp_1.pdf','-transparent');

