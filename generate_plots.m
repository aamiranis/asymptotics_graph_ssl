addpath(genpath('~/matlab/plotting'));

%% Density separable

clear;
load results/sep_exp_1.mat;

figure1 = figure;

font_size = 16;

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


%% Density nonseparable

clear;

load results/nonsep_exp_1.mat;

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

zlim([0 0.6]);
set(gca,'ZTick',0:0.2:0.6);
set(gca,'YTick',-2:2:2);
set(gca,'XTick',-3:1:3);
set(gca, 'FontSize', font_size);
xlabel('x','FontSize', font_size);
ylabel('y','FontSize', font_size);
zlabel('p(x,y)','FontSize', font_size);

export_fig('plots/nonsep_density.pdf','-transparent');

%% Boundaries: separable

clear;
load results/sep_exp_1.mat;

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

%% Bandwidth convergence

clear;
load results/sep_exp_1.mat;

emp_means = mean_bw_list;
emp_std_devs = std_bw_list;

theoretical = sup_cut';

load results/nonsep_exp_1.mat;

emp_means = [emp_means; mean_bw_list;];
emp_std_devs = [emp_std_devs; std_bw_list;];

theoretical = [theoretical; sup_cut;];

font_size = 16;

figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.45]);

ax1 = bar([emp_means + emp_std_devs/2 theoretical],'LineWidth',1);
set(ax1(1),'FaceColor',[0.6 0.6 0.6]);
set(ax1(2),'FaceColor',[1 0.5 0.5]);
ylim([0 1]);
set(gca,'FontSize',font_size);
set(gca,'XTickLabel',{'','','','','',''});
annotation(figure1,'textbox',[0.22 0.055 0.018 0.048],'String',{'\partialS_1'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.33 0.055 0.018 0.048],'String',{'\partialS_2'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.44 0.055 0.018 0.048],'String',{'\partialS_3'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.55 0.055 0.018 0.048],'String',{'\partialS_4'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.66 0.055 0.018 0.048],'String',{'\partialS_5'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
annotation(figure1,'textbox',[0.77 0.075 0.018 0.048],'String',{'\partialA'},...
    'FontSize',font_size,'FitBoxToText','off','EdgeColor','none');
set(gca,'YTick',0:0.2:1);
hold all;

ax2 = bar([emp_means - emp_std_devs/2 theoretical],'LineWidth',1);
set(ax2(1),'FaceColor',[0.9 0.9 0.9]);
set(ax2(2),'FaceColor',[1 0.5 0.5]);

ax3 = bar([emp_means theoretical],'LineWidth',1);
set(ax3(1),'FaceColor','none');
set(ax3(2),'FaceColor',[1 0.5 0.5]);

ylabel('Bandwidth');

legend1 = legend([ax2(1) ax2(2)],{'Empirical', 'Theoretical'});
set(legend1,...
    'Position',[0.15 0.77 0.28 0.08]);

export_fig('plots/exp_1.pdf','-transparent');

%% Bandwidth convergence with n

clear;
load results/sep_exp_2.mat;

std_devs = std_bw_list(1:end-1,:);

load results/nonsep_exp_2.mat;
std_devs = [std_devs std_bw_list(1:end-1)'];

font_size = 16;

figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.3]);

plot(N_list(1:end-1), std_devs, '--o', 'LineWidth',2);
xlabel('n','FontSize',font_size);
ylabel('std. dev.','FontSize',font_size);
set(gca,'FontSize',font_size);
xlim([0 3500]);
set(gca,'XTick',0:500:3000);
ylim([0 0.1]);

legend1 = legend('\partial S_1','\partial S_2','\partial S_3','\partial S_4','\partial S_5','\partial A', 'Location', 'northeastoutside');
set(legend1, 'FontSize', 14);

axis_position = get(gca, 'Position');
position = get(legend1, 'Position');
position(2) = position(2) - 0.09;
position(3) = position(3) + 0.025;
set(legend1, 'Position', position);
set(gca, 'Position', axis_position);

export_fig('plots/exp_2.pdf','-transparent');

%% Reconstruction label complexity: separable

clear;
load results/sep_exp_3.mat;

font_size = 16;

figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.3]);

plot(label_percent_list, mean_error_list_3, 'k-^', [lc_cut(2) lc_cut(2)], [0 0.1], 'r--','LineWidth',2);
ylabel('Mean error','FontSize',font_size);
xlabel('Fraction of labeled examples','FontSize',font_size);
set(gca,'FontSize',font_size);
ylim([0 0.1]);
set(gca,'YTick',0:0.02:0.1);
export_fig('plots/exp_3_sep.pdf','-transparent');


%% Reconstruction label complexity: nonseparable

clear;
load results/nonsep_exp_3.mat;

font_size = 16;

figure1 = figure;
scr = get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(figure1, 'Position', [scr(3)*0.25 scr(4)*0.3 scr(3)*0.5 scr(4)*0.3]);

plot(label_percent_list, mean_error_list, 'k-^', [lc_cut lc_cut], [0 0.1], 'r--','LineWidth',2);
ylabel('Mean error','FontSize',font_size);
xlabel('Fraction of labeled examples','FontSize',font_size);
set(gca,'FontSize',font_size);
ylim([0 0.05]);
set(gca,'YTick',0:0.01:0.5);
export_fig('plots/exp_3_nonsep.pdf','-transparent');
