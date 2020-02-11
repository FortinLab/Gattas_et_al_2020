% steps: run A006 for all animals then run this.
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
animal_length =6;
clear cond1 cond2 cond3 cond4 bar_plot_4cond_values_per_anim

% load group data
if strcmp('novel1',task)||strcmp('novel2',task)
    load(['bar_plot_4cond_values_per_anim_' band '_nonmatched_' task])
elseif strcmp('welltrained',task)
    if strcmp( 'yes', match)
        load(['bar_plot_4cond_values_per_anim_' band '_matched'])
    else
        load(['bar_plot_4cond_values_per_anim_' band '_' band_value '_nonmatched_' task])
    end
end


if strcmp('novel1', task) || strcmp('novel2', task)
    %inseq+
    cond1 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
        bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)]; % (elec, anim)
        
        
        %outseq-
        cond2 =  [bar_plot_4cond_values_per_anim{1}(:,2) bar_plot_4cond_values_per_anim{2}(:,2) ...
        bar_plot_4cond_values_per_anim{3}(:,2)  bar_plot_4cond_values_per_anim{4}(:,2)]; % (elec, anim)
        
        
        %inseq-
        cond3= [bar_plot_4cond_values_per_anim{1}(:,3)  bar_plot_4cond_values_per_anim{2}(:,3) ...
        bar_plot_4cond_values_per_anim{3}(:,3)  bar_plot_4cond_values_per_anim{4}(:,3)]; % (elec, anim)
        
        
        %outseq-
        cond4= [ bar_plot_4cond_values_per_anim{1}(:,4) bar_plot_4cond_values_per_anim{2}(:,4) ...
        bar_plot_4cond_values_per_anim{3}(:,4)  bar_plot_4cond_values_per_anim{4}(:,4)]; % (elec, anim)
        
else
    if animal_length ==6
        %inseq+
        cond1 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
            bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)...
            bar_plot_4cond_values_per_anim{5}(:,1) bar_plot_4cond_values_per_anim{6}(:,1)]; % (elec, anim)
        
        
        %outseq-
        cond2 =  [bar_plot_4cond_values_per_anim{1}(:,2) bar_plot_4cond_values_per_anim{2}(:,2) ...
            bar_plot_4cond_values_per_anim{3}(:,2)  bar_plot_4cond_values_per_anim{4}(:,2)...
            bar_plot_4cond_values_per_anim{5}(:,2) bar_plot_4cond_values_per_anim{6}(:,2)]; % (elec, anim)
        
        
        %inseq-
        cond3= [bar_plot_4cond_values_per_anim{1}(:,3)  bar_plot_4cond_values_per_anim{2}(:,3) ...
            bar_plot_4cond_values_per_anim{3}(:,3)  bar_plot_4cond_values_per_anim{4}(:,3)...
            bar_plot_4cond_values_per_anim{5}(:,3) bar_plot_4cond_values_per_anim{6}(:,3)]; % (elec, anim)
        
        
        %outseq-
        cond4= [ bar_plot_4cond_values_per_anim{1}(:,4) bar_plot_4cond_values_per_anim{2}(:,4) ...
            bar_plot_4cond_values_per_anim{3}(:,4)  bar_plot_4cond_values_per_anim{4}(:,4)...
            bar_plot_4cond_values_per_anim{5}(:,4) bar_plot_4cond_values_per_anim{6}(:,4)]; % (elec, anim)
    elseif animal_length ==5
        %inseq+
        cond1 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
            bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)...
            bar_plot_4cond_values_per_anim{5}(:,1) ]; % (elec, anim)
        
        
        %outseq-
        cond2 =  [bar_plot_4cond_values_per_anim{1}(:,2) bar_plot_4cond_values_per_anim{2}(:,2) ...
            bar_plot_4cond_values_per_anim{3}(:,2)  bar_plot_4cond_values_per_anim{4}(:,2)...
            bar_plot_4cond_values_per_anim{5}(:,2) ]; % (elec, anim)
        
        
        %inseq-
        cond3= [bar_plot_4cond_values_per_anim{1}(:,3)  bar_plot_4cond_values_per_anim{2}(:,3) ...
            bar_plot_4cond_values_per_anim{3}(:,3)  bar_plot_4cond_values_per_anim{4}(:,3)...
            bar_plot_4cond_values_per_anim{5}(:,3)]; % (elec, anim)
        
        
        %outseq-
        cond4= [ bar_plot_4cond_values_per_anim{1}(:,4) bar_plot_4cond_values_per_anim{2}(:,4) ...
            bar_plot_4cond_values_per_anim{3}(:,4)  bar_plot_4cond_values_per_anim{4}(:,4)...
            bar_plot_4cond_values_per_anim{5}(:,4) ]; % (elec, anim)
    end
end
%%
mn = -.6
mx = .6
figure

% elec 1
subplot(1, 4, 1)
bar_vectora1 = [nanmean(cond1(1,:)) nanmean(cond4(1,:)) nanmean(cond2(1,:)) nanmean(cond3(1,:))];
bar(bar_vectora1)
hold on
errorbar(1:4,bar_vectora1,[nanstd(cond1(1,:,1))/sqrt(animal_length) nanstd(cond4(1,:,1))/sqrt(animal_length) nanstd(cond2(1,:,1))/sqrt(animal_length) nanstd(cond3(1,:,1))/sqrt(animal_length)], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 2
subplot(1, 4, 2)
bar_vectora2 = [nanmean(cond1(2,:)) nanmean(cond4(2,:)) nanmean(cond2(2,:)) nanmean(cond3(2,:))];
bar(bar_vectora2)
hold on
errorbar(1:4,bar_vectora2,[nanstd(cond1(2,:,1))/sqrt(animal_length) nanstd(cond4(2,:,1))/sqrt(animal_length) nanstd(cond2(2,:,1))/sqrt(animal_length) nanstd(cond3(2,:,1))/sqrt(animal_length)], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 3
subplot(1, 4, 3)
bar_vectora3 = [nanmean(cond1(3,:)) nanmean(cond4(3,:)) nanmean(cond2(3,:)) nanmean(cond3(3,:))];
bar(bar_vectora3)
hold on
errorbar(1:4,bar_vectora3,[nanstd(cond1(3,:,1))/sqrt(animal_length) nanstd(cond4(3,:,1))/sqrt(animal_length) nanstd(cond2(3,:,1))/sqrt(animal_length) nanstd(cond3(3,:,1))/sqrt(animal_length)], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 4
subplot(1, 4, 4)
bar_vectora4 = [nanmean(cond1(4,:)) nanmean(cond4(4,:)) nanmean(cond2(4,:)) nanmean(cond3(4,:))];
bar(bar_vectora4)
hold on
errorbar(1:4,bar_vectora4,[nanstd(cond1(4,:,1))/sqrt(animal_length) nanstd(cond4(4,:,1))/sqrt(animal_length) nanstd(cond2(4,:,1))/sqrt(animal_length) nanstd(cond3(4,:,1))/sqrt(animal_length)], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'InSeq +' 'OutSeq -' 'OutSeq +' 'InSeq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
supertitle([band ' band: ' band_value], 'Interpreter', 'none')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

%%
figure
bar_vector = [(nanmean(cond1(1,:))- nanmean(cond4(1,:))) (nanmean(cond1(2,:))- nanmean(cond4(2,:))) (nanmean(cond1(3,:))- nanmean(cond4(3,:))) (nanmean(cond1(4,:))- nanmean(cond4(4,:)))];
bar(bar_vector)

figure
bar_vector = [(nanmean(cond1(1,:))- nanmean(cond2(1,:))) (nanmean(cond1(2,:))- nanmean(cond2(2,:))) (nanmean(cond1(3,:))- nanmean(cond2(3,:))) (nanmean(cond1(4,:))- nanmean(cond2(4,:)))];
bar(bar_vector)


figure
bar_vector = [(nanmean(cond4(1,:))- nanmean(cond2(1,:))) (nanmean(cond4(2,:))- nanmean(cond2(2,:))) (nanmean(cond4(3,:))- nanmean(cond2(3,:))) (nanmean(cond4(4,:))- nanmean(cond2(4,:)))];
bar(bar_vector)

figure
bar_vector = [(nanmean(cond4(1,:))- nanmean(cond3(1,:))) (nanmean(cond4(2,:))- nanmean(cond3(2,:))) (nanmean(cond4(3,:))- nanmean(cond3(3,:))) (nanmean(cond4(4,:))- nanmean(cond3(4,:)))];
bar(bar_vector)

%% another freq to do a ratio
clear cond1 cond2 cond3 cond4 bar_plot_4cond_values_per_anim

band = 'theta'
band_value = '9-12'
load(['bar_plot_4cond_values_per_anim_' band '_' band_value '_nonmatched_' task])


if strcmp('novel1', task) || strcmp('novel2', task)
    %inseq+
    cond1 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
        bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)]; % (elec, anim)
        
        
        %outseq-
        cond2 =  [bar_plot_4cond_values_per_anim{1}(:,2) bar_plot_4cond_values_per_anim{2}(:,2) ...
        bar_plot_4cond_values_per_anim{3}(:,2)  bar_plot_4cond_values_per_anim{4}(:,2)]; % (elec, anim)
        
        
        %inseq-
        cond3= [bar_plot_4cond_values_per_anim{1}(:,3)  bar_plot_4cond_values_per_anim{2}(:,3) ...
        bar_plot_4cond_values_per_anim{3}(:,3)  bar_plot_4cond_values_per_anim{4}(:,3)]; % (elec, anim)
        
        
        %outseq-
        cond4= [ bar_plot_4cond_values_per_anim{1}(:,4) bar_plot_4cond_values_per_anim{2}(:,4) ...
        bar_plot_4cond_values_per_anim{3}(:,4)  bar_plot_4cond_values_per_anim{4}(:,4)]; % (elec, anim)
        
else
%inseq+
cond1 = [ bar_plot_4cond_values_per_anim{1}(:,1) bar_plot_4cond_values_per_anim{2}(:,1) ...
    bar_plot_4cond_values_per_anim{3}(:,1)  bar_plot_4cond_values_per_anim{4}(:,1)...
    bar_plot_4cond_values_per_anim{5}(:,1) bar_plot_4cond_values_per_anim{6}(:,1)]; % (elec, anim)


%outseq-
cond2 =  [bar_plot_4cond_values_per_anim{1}(:,2) bar_plot_4cond_values_per_anim{2}(:,2) ...
    bar_plot_4cond_values_per_anim{3}(:,2)  bar_plot_4cond_values_per_anim{4}(:,2)...
    bar_plot_4cond_values_per_anim{5}(:,2) bar_plot_4cond_values_per_anim{6}(:,2)]; % (elec, anim)


%inseq-
cond3= [bar_plot_4cond_values_per_anim{1}(:,3)  bar_plot_4cond_values_per_anim{2}(:,3) ...
    bar_plot_4cond_values_per_anim{3}(:,3)  bar_plot_4cond_values_per_anim{4}(:,3)...
    bar_plot_4cond_values_per_anim{5}(:,3) bar_plot_4cond_values_per_anim{6}(:,3) ]; % (elec, anim)


%outseq-
cond4= [ bar_plot_4cond_values_per_anim{1}(:,4) bar_plot_4cond_values_per_anim{2}(:,4) ...
    bar_plot_4cond_values_per_anim{3}(:,4)  bar_plot_4cond_values_per_anim{4}(:,4)...
    bar_plot_4cond_values_per_anim{5}(:,4) bar_plot_4cond_values_per_anim{6}(:,4)]; % (elec, anim)
end

%%
mn = -10;
mx = 90;
figure

% elec 1
subplot(1, 4, 1)
bar_vectorb = [nanmean(cond1(1,:)) nanmean(cond4(1,:)) nanmean(cond2(1,:)) nanmean(cond3(1,:))];
bar(bar_vectora1./bar_vectorb)
hold on
%errorbar(1:4,bar_vectora./bar_vectorb,[nanstd(cond1(1,:,1)) nanstd(cond4(1,:,1)) nanstd(cond2(1,:,1)) nanstd(cond3(1,:,1))], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 2
subplot(1, 4, 2)
bar_vectorb = [nanmean(cond1(2,:)) nanmean(cond4(2,:)) nanmean(cond2(2,:)) nanmean(cond3(2,:))];
bar(bar_vectora2./bar_vectorb)
hold on
%errorbar(1:4,bar_vector,[nanstd(cond1(2,:,1)) nanstd(cond4(2,:,1)) nanstd(cond2(2,:,1)) nanstd(cond3(2,:,1))], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 3
subplot(1, 4, 3)
bar_vectorb = [nanmean(cond1(3,:)) nanmean(cond4(3,:)) nanmean(cond2(3,:)) nanmean(cond3(3,:))];
bar(bar_vectora3./bar_vectorb)
hold on
%errorbar(1:4,bar_vector,[nanstd(cond1(3,:,1)) nanstd(cond4(3,:,1)) nanstd(cond2(3,:,1)) nanstd(cond3(3,:,1))], 'rx')
ylim([mn mx])

set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 4
subplot(1, 4, 4)
bar_vectorb = [nanmean(cond1(4,:)) nanmean(cond4(4,:)) nanmean(cond2(4,:)) nanmean(cond3(4,:))];
bar(bar_vectora4./bar_vectorb)
hold on
%errorbar(1:4,bar_vector,[nanstd(cond1(4,:,1)) nanstd(cond4(4,:,1)) nanstd(cond2(4,:,1)) nanstd(cond3(4,:,1))], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
supertitle(['beta:theta ratio'], 'Interpreter', 'none')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

%% plot only inseq corr vs. outseq incorr
mn = -.45;
mx = .2;
figure

% elec 1
subplot(1, 4, 1)
bar_vectora1 = [nanmean(cond1(1,:)) nanmean(cond4(1,:)) ];
bar(bar_vectora1)
hold on
errorbar(1:2,bar_vectora1,[nanstd(cond1(1,:,1))/sqrt(animal_length) nanstd(cond4(1,:,1))/sqrt(animal_length)], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 2
subplot(1, 4, 2)
bar_vectora2 = [nanmean(cond1(2,:)) nanmean(cond4(2,:))];
bar(bar_vectora2)
hold on
errorbar(1:2,bar_vectora2,[nanstd(cond1(2,:,1))/sqrt(animal_length) nanstd(cond4(2,:,1))/sqrt(animal_length) ], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:2, 'XTickLabel', {'inseq +' 'otseq -' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 3
subplot(1, 4, 3)
bar_vectora3 = [nanmean(cond1(3,:)) nanmean(cond4(3,:))];
bar(bar_vectora3)
hold on
errorbar(1:2,bar_vectora3,[nanstd(cond1(3,:,1))/sqrt(animal_length) nanstd(cond4(3,:,1))/sqrt(animal_length)  ], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

% elec 4
subplot(1, 4, 4)
bar_vectora4 = [nanmean(cond1(4,:)) nanmean(cond4(4,:))];
bar(bar_vectora4)
hold on
errorbar(1:2,bar_vectora4,[nanstd(cond1(4,:,1))/sqrt(animal_length) nanstd(cond4(4,:,1))/sqrt(animal_length) ], 'rx')
ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' }, 'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
