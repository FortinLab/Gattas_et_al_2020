band = 'beta'
animal_length = 4
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')

%welltrained
load(['bar_plot_4cond_values_per_anim_' band '_nonmatched'])
cond1 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
    bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)]; % (elec, anim)

load(['bar_plot_4cond_values_per_anim_' band '_nonmatched_novel2'])
%novel 1
cond2 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
    bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)]; % (elec, anim)

load(['bar_plot_4cond_values_per_anim_' band '_nonmatched_novel1'])

%novel 2
cond3 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
    bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)]; % (elec, anim)

%%

mn = -.2
mx = .6
figure
% elec 1
subplot(1, 4, 1)
bar_vector = [nanmean(cond1(1,:)) nanmean(cond2(1,:)) nanmean(cond3(1,:))];
bar(bar_vector)
hold on
errorbar(1:3,bar_vector,[nanstd(cond1(1,:,1))/sqrt(animal_length) nanstd(cond2(1,:,1))/sqrt(animal_length) nanstd(cond3(1,:,1))/sqrt(animal_length)], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:3, 'XTickLabel', {'welltrained' 'novel 2' 'novel 1'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 2
subplot(1, 4, 2)
bar_vector = [nanmean(cond1(2,:)) nanmean(cond2(2,:)) nanmean(cond3(2,:))];
bar(bar_vector)
hold on
errorbar(1:3,bar_vector,[nanstd(cond1(2,:,1))/sqrt(animal_length) nanstd(cond2(2,:,1))/sqrt(animal_length) nanstd(cond3(2,:,1))/sqrt(animal_length)], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:3, 'XTickLabel', {'welltrained' 'novel 2' 'novel 1'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 3
subplot(1, 4, 3)
bar_vector = [nanmean(cond1(3,:)) nanmean(cond2(3,:)) nanmean(cond3(3,:)) ];
bar(bar_vector)
hold on
errorbar(1:3,bar_vector,[nanstd(cond1(3,:,1))/sqrt(animal_length) nanstd(cond2(3,:,1))/sqrt(animal_length) nanstd(cond3(3,:,1))/sqrt(animal_length) ], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:3, 'XTickLabel', {'welltrained' 'novel 2' 'novel 1'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 4
subplot(1, 4, 4)
bar_vector = [nanmean(cond1(4,:)) nanmean(cond2(4,:)) nanmean(cond3(4,:))];
bar(bar_vector)
hold on
errorbar(1:3,bar_vector,[nanstd(cond1(4,:,1))/sqrt(animal_length) nanstd(cond2(4,:,1))/sqrt(animal_length) nanstd(cond3(4,:,1))/sqrt(animal_length) ], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:3, 'XTickLabel', {'welltrained' 'novel 2' 'novel 1'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
