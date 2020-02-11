% first run A006 stim locked bar plots matched trial num acrs chan anim

%% start and end time idx
cd('D:\Gattas\ephys_data_final\welltrained\group_plots')
if strcmp('p_o',lock)
    trial_num = 51;
else
    trial_num = 49;
end
cond_1_all_animals = nan(size(cond_4,1), size(cond_4,2), trial_num , length(direc_chan));
cond_2_all_animals = nan(size(cond_4,1), size(cond_4,2), trial_num , length(direc_chan));
cond_3_all_animals = nan(size(cond_4,1), size(cond_4,2), trial_num , length(direc_chan));
cond_4_all_animals = nan(size(cond_4,1), size(cond_4,2), trial_num , length(direc_chan));

if strcmp('p_o', lock)
cond1 = [9 9 9 8 8 8 ]
cond2 = [8 8 8 9 9 9 ]
cond3 = [0 17 1 6 17 10]
cond4 = [3 5 8 15 14 6]
else
cond1 = [8 9 8 8 8 8]
cond2 = [9 8 8 8 8 8]
cond3 = [0 15 1 7 16 10]
cond4 = [2 4 9 15 13 6] 
end

anim=1;
if strcmp('p_o', lock)
    load(['anim' num2str(anim) '_exact_match_' lock])
else
    load(['anim' num2str(anim) '_exact_match_' lock])
end

for chan = 1:4
cond_1_all_animals(:,:,1:cond1(anim), chan) = cond_1(:,:,:,chan);
cond_2_all_animals(:,:,1:cond2(anim), chan) = cond_2(:,:,:,chan);
cond_4_all_animals(:,:,1:cond4(anim), chan) = cond_4(:,:,:,chan);
end

for anim = 2:6
    if strcmp('p_o', lock)
        load(['anim' num2str(anim) '_exact_match' ])
    else
        load(['anim' num2str(anim) '_exact_match_' lock])
    end
    
    for chan = 1:4
        cond_1_all_animals(:,:,(sum(cond1(1:anim-1))+1):(sum(cond1(1:anim-1)))+cond1(anim),chan) = cond_1(:,:,:,chan);
        cond_2_all_animals(:,:,(sum(cond2(1:anim-1))+1):(sum(cond2(1:anim-1)))+cond2(anim),chan) = cond_2(:,:,:,chan);
        cond_3_all_animals(:,:,(sum(cond3(1:anim-1))+1):(sum(cond3(1:anim-1)))+cond3(anim),chan) = cond_3(:,:,:,chan);
        cond_4_all_animals(:,:,(sum(cond4(1:anim-1))+1):(sum(cond4(1:anim-1)))+cond4(anim),chan) = cond_4(:,:,:,chan);
    end
end


%% index time and freq of interest, then average
% freq band
freq_range= freq>19 & freq<36;  % experimental cond
band_value = '19-36 hz';
min_times = [ 483 178  527  284 27 318] % all trials
min_times = [481] % outseq corr
start_time = 0.5*fs - min(min_times);

cond_1_mean_all_trials = []
cond_2_mean_all_trials = []
cond_3_mean_all_trials = []
cond_4_mean_all_trials = []

%% correct for std then bar plot
start_time = 250; % 250ms
end_time   = 500;   % 600ms
cond_1_mean_each_trial = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);
cond_1_std_each_trial  = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);
cond_2_mean_each_trial = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);
cond_2_std_each_trial  = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);
cond_3_mean_each_trial = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);
cond_3_std_each_trial  = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);
cond_4_mean_each_trial = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);
cond_4_std_each_trial  = nan(size(cond_1_all_animals,3),length(start_time:end_time),4);

for chan = 1:4 % loop thru chans
   for trial = 1:size(cond_1_all_animals,3) % loop thru trials
       % cond 1
       cond_1_mean_each_trial(trial,:,chan) = mean(cond_1_all_animals(freq_range,start_time:end_time,trial,chan),1);
       cond_1_std_each_trial(trial,:,chan) = std(cond_1_all_animals(freq_range,start_time:end_time,trial,chan),0,1);
       % cond 2
       cond_2_mean_each_trial(trial,:,chan) = mean(cond_2_all_animals(freq_range,start_time:end_time,trial,chan),1);
       cond_2_std_each_trial(trial,:,chan) = std(cond_2_all_animals(freq_range,start_time:end_time,trial,chan),0,1);
       % cond 3
       cond_3_mean_each_trial(trial,:,chan) = mean(cond_3_all_animals(freq_range,start_time:end_time,trial,chan),1);
       cond_3_std_each_trial(trial,:,chan) = std(cond_3_all_animals(freq_range,start_time:end_time,trial,chan),0,1);
       % cond 4
       cond_4_mean_each_trial(trial,:,chan) = mean(cond_4_all_animals(freq_range,start_time:end_time,trial,chan),1);
       cond_4_std_each_trial(trial,:,chan) = std(cond_4_all_animals(freq_range,start_time:end_time,trial,chan),0,1);
   end
end
% divide the mean by std
cond_1_corrected = squeeze(mean(cond_1_mean_each_trial,1))'./squeeze(mean(cond_1_std_each_trial,1))';
cond_2_corrected = squeeze(mean(cond_2_mean_each_trial,1))'./squeeze(mean(cond_2_std_each_trial,1))';
cond_3_corrected = squeeze(mean(cond_3_mean_each_trial,1))'./squeeze(mean(cond_3_std_each_trial,1))';
cond_4_corrected = squeeze(mean(cond_4_mean_each_trial,1))'./squeeze(mean(cond_4_std_each_trial,1))';

cond_1_noncorrected = squeeze(mean(cond_1_mean_each_trial,1))';
cond_2_noncorrected = squeeze(mean(cond_2_mean_each_trial,1))';
cond_3_noncorrected = squeeze(mean(cond_3_mean_each_trial,1))';
cond_4_noncorrected = squeeze(mean(cond_4_mean_each_trial,1))';



%% time traces with corrected means
figure;
elec = 1
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_corrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_corrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_corrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_corrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

elec = 2
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_corrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_corrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_corrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_corrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

elec = 3
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_corrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_corrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_corrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_corrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

elec = 4
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_corrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_corrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_corrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_corrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')


%% time traces with non-corrected means
figure;
elec = 1
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_noncorrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_noncorrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_noncorrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_noncorrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

elec = 2
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_noncorrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_noncorrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_noncorrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_noncorrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

elec = 3
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_noncorrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_noncorrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_noncorrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_noncorrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

elec = 4
subplot(1,4,elec)
hold on
ylim ([-1.5 1.5])
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_1_noncorrected(elec,:), 'r', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_4_noncorrected(elec,:), 'm', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_2_noncorrected(elec,:), 'b', 'LineWidth', 1.8)
plot(linspace(eventWindow(1), eventWindow(2),size(cond_1_mean_each_trial,2)), cond_3_noncorrected(elec,:), 'g', 'LineWidth', 1.8)
plot(zeros(1,length(1:1000)), linspace(-1.5,1.5, length(1:1000)),'Color','k','LineWidth',1)
legend({'inseq +', 'otseq -', 'otseq +', 'inseq -'},'FontSize',14, 'FontWeight', 'bold')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

%% bar plots with corrected means
figure
mn = -.3
mx = 1
% elec 1
elec = 1
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_corrected(elec,:)) mean(cond_4_corrected(elec,:)) mean(cond_2_corrected(elec,:)) mean(cond_3_corrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_corrected(elec,:)) std(cond_4_corrected(elec,:)) std(cond_2_corrected(elec,:)) std(cond_3_corrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_corrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])



elec =2
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_corrected(elec,:)) mean(cond_4_corrected(elec,:)) mean(cond_2_corrected(elec,:)) mean(cond_3_corrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_corrected(elec,:)) std(cond_4_corrected(elec,:)) std(cond_2_corrected(elec,:)) std(cond_3_corrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_corrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

elec = 3
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_corrected(elec,:)) mean(cond_4_corrected(elec,:)) mean(cond_2_corrected(elec,:)) mean(cond_3_corrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_corrected(elec,:)) std(cond_4_corrected(elec,:)) std(cond_2_corrected(elec,:)) std(cond_3_corrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_corrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

elec = 4
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_corrected(elec,:)) mean(cond_4_corrected(elec,:)) mean(cond_2_corrected(elec,:)) mean(cond_3_corrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_corrected(elec,:)) std(cond_4_corrected(elec,:)) std(cond_2_corrected(elec,:)) std(cond_3_corrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_corrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

%% bar plots with noncorrected means
figure
mn = -.3
mx = 1
% elec 1
elec = 1
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_noncorrected(elec,:)) mean(cond_4_noncorrected(elec,:)) mean(cond_2_noncorrected(elec,:)) mean(cond_3_noncorrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_noncorrected(elec,:)) std(cond_4_noncorrected(elec,:)) std(cond_2_noncorrected(elec,:)) std(cond_3_noncorrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_noncorrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45,'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

elec =2
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_noncorrected(elec,:)) mean(cond_4_noncorrected(elec,:)) mean(cond_2_noncorrected(elec,:)) mean(cond_3_noncorrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_noncorrected(elec,:)) std(cond_4_noncorrected(elec,:)) std(cond_2_noncorrected(elec,:)) std(cond_3_noncorrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_noncorrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45,'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

elec = 3
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_noncorrected(elec,:)) mean(cond_4_noncorrected(elec,:)) mean(cond_2_noncorrected(elec,:)) mean(cond_3_noncorrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_noncorrected(elec,:)) std(cond_4_noncorrected(elec,:)) std(cond_2_noncorrected(elec,:)) std(cond_3_noncorrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_noncorrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45,'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

elec = 4
subplot(1, 4, elec)
bar_vectora1 = [mean(cond_1_noncorrected(elec,:)) mean(cond_4_noncorrected(elec,:)) mean(cond_2_noncorrected(elec,:)) mean(cond_3_noncorrected(elec,:)) ];
bar(bar_vectora1)
hold on
bar_vectora1std = [std(cond_1_noncorrected(elec,:)) std(cond_4_noncorrected(elec,:)) std(cond_2_noncorrected(elec,:)) std(cond_3_noncorrected(elec,:))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_noncorrected,2), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45,'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

%% dont correct for std
for chan = 1:4
    cond_1_mean_all_trials(:,chan)=  mean(squeeze(mean(cond_1_all_animals(freq_range,start_time:end_time,:,chan),1)),1);
    cond_2_mean_all_trials(:,chan)=  mean(squeeze(mean(cond_2_all_animals(freq_range,start_time:end_time,:,chan),1)),1);
    cond_3_mean_all_trials(:,chan)=  mean(squeeze(mean(cond_3_all_animals(freq_range,start_time:end_time,:,chan),1)),1);
    cond_4_mean_all_trials(:,chan)=  mean(squeeze(mean(cond_4_all_animals(freq_range,start_time:end_time,:,chan),1)),1);
end

cond_1_mean_trace = []
cond_2_mean_trace = []
cond_3_mean_trace = []
cond_4_mean_trace = []
for chan = 1:4
    cond_1_mean_trace(:,:,chan)= squeeze(mean(cond_1_all_animals(freq_range,start_time:end_time,:,chan),1))';
    cond_2_mean_trace(:,:,chan)=  squeeze(mean(cond_2_all_animals(freq_range,start_time:end_time,:,chan),1))';
    cond_3_mean_trace(:,:,chan)=  squeeze(mean(cond_3_all_animals(freq_range,start_time:end_time,:,chan),1))';
    cond_4_mean_trace(:,:,chan)= squeeze(mean(cond_4_all_animals(freq_range,start_time:end_time,:,chan),1))';
end


%% generate bar plots

figure
mn = -.5
mx = 1
% elec 1
elec = 1
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_1_mean_all_trials(:,elec)) nanmean(cond_4_mean_all_trials(:,elec)) nanmean(cond_2_mean_all_trials(:,elec)) nanmean(cond_3_mean_all_trials(:,elec))];
bar(bar_vectora1)
hold on
bar_vectora1std = [nanstd(cond_1_mean_all_trials(:,elec)) nanstd(cond_4_mean_all_trials(:,elec)) nanstd(cond_2_mean_all_trials(:,elec)) nanstd(cond_3_mean_all_trials(:,elec))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_mean_all_trials,1), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

% elec 2
elec = 2
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_1_mean_all_trials(:,elec)) nanmean(cond_4_mean_all_trials(:,elec)) nanmean(cond_2_mean_all_trials(:,elec)) nanmean(cond_3_mean_all_trials(:,elec))];
bar(bar_vectora1)
hold on
bar_vectora1std = [nanstd(cond_1_mean_all_trials(:,elec)) nanstd(cond_4_mean_all_trials(:,elec)) nanstd(cond_2_mean_all_trials(:,elec)) nanstd(cond_3_mean_all_trials(:,elec))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_mean_all_trials,1), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

% elec 3
elec = 3
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_1_mean_all_trials(:,elec)) nanmean(cond_4_mean_all_trials(:,elec)) nanmean(cond_2_mean_all_trials(:,elec)) nanmean(cond_3_mean_all_trials(:,elec))];
bar(bar_vectora1)
hold on
bar_vectora1std = [nanstd(cond_1_mean_all_trials(:,elec)) nanstd(cond_4_mean_all_trials(:,elec)) nanstd(cond_2_mean_all_trials(:,elec)) nanstd(cond_3_mean_all_trials(:,elec))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_mean_all_trials,1), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

% elec 4
elec = 4
subplot(1, 4, elec)
bar_vectora1 = [nanmean(cond_1_mean_all_trials(:,elec)) nanmean(cond_4_mean_all_trials(:,elec)) nanmean(cond_2_mean_all_trials(:,elec)) nanmean(cond_3_mean_all_trials(:,elec))];
bar(bar_vectora1)
hold on
bar_vectora1std = [nanstd(cond_1_mean_all_trials(:,elec)) nanstd(cond_4_mean_all_trials(:,elec)) nanstd(cond_2_mean_all_trials(:,elec)) nanstd(cond_3_mean_all_trials(:,elec))];
errorbar(1:4,bar_vectora1,bar_vectora1std/size(cond_1_mean_all_trials,1), 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'otseq -' 'otseq +' 'inseq -'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])

%% time domain plots
win = .1*fs %100 ms smoothing window

figure
mx = .9
mn = -.9
elec = 1
begin_idx = -.5
end_idx = .5
hold on
subplot(4, 1, elec)
stdshade(cond_1_mean_trace(:,:,elec),.2,'r',linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ),[])
stdshade(cond_4_mean_trace(:,:,elec),.2,'m',linspace(begin_idx, end_idx,size(cond_4_mean_trace,2) ),[])
stdshade(cond_2_mean_trace(:,:,elec),.2,'b',linspace(begin_idx, end_idx,size(cond_2_mean_trace,2) ),[])
stdshade(cond_3_mean_trace(:,:,elec),.2,'g',linspace(begin_idx, end_idx,size(cond_3_mean_trace,2) ),[])
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
legend({'inseq +', 'SEM', 'otseq -', 'SEM', 'otseq +', 'SEM', 'inseq -', 'SEM'},'FontSize',10, 'FontWeight', 'bold')
ylim ([mn mx])



elec = 2
subplot(4, 1, elec)

hold on
stdshade(cond_1_mean_trace(:,:,elec),.2,'r',linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ),[])
stdshade(cond_4_mean_trace(:,:,elec),.2,'m',linspace(begin_idx, end_idx,size(cond_4_mean_trace,2) ),[])
stdshade(cond_2_mean_trace(:,:,elec),.2,'b',linspace(begin_idx, end_idx,size(cond_2_mean_trace,2) ),[])
stdshade(cond_3_mean_trace(:,:,elec),.2,'g',linspace(begin_idx, end_idx,size(cond_3_mean_trace,2) ),[])
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
%legend({'inseq +', 'SEM', 'otseq -', 'SEM', 'otseq +', 'SEM', 'inseq -', 'SEM'},'FontSize',10, 'FontWeight', 'bold')
ylim ([mn mx])


elec = 3
subplot(4, 1, elec)

hold on
stdshade(cond_1_mean_trace(:,:,elec),.2,'r',linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ),[])
stdshade(cond_4_mean_trace(:,:,elec),.2,'m',linspace(begin_idx, end_idx,size(cond_4_mean_trace,2) ),[])
stdshade(cond_2_mean_trace(:,:,elec),.2,'b',linspace(begin_idx, end_idx,size(cond_2_mean_trace,2) ),[])
stdshade(cond_3_mean_trace(:,:,elec),.2,'g',linspace(begin_idx, end_idx,size(cond_3_mean_trace,2) ),[])
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
%legend({'inseq +', 'SEM', 'otseq -', 'SEM', 'otseq +', 'SEM', 'inseq -', 'SEM'},'FontSize',10, 'FontWeight', 'bold')
ylim ([mn mx])


elec = 4
subplot(4, 1, elec)

hold on
stdshade(cond_1_mean_trace(:,:,elec),.2,'r',linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ),[])
stdshade(cond_4_mean_trace(:,:,elec),.2,'m',linspace(begin_idx, end_idx,size(cond_4_mean_trace,2) ),[])
stdshade(cond_2_mean_trace(:,:,elec),.2,'b',linspace(begin_idx, end_idx,size(cond_2_mean_trace,2) ),[])
stdshade(cond_3_mean_trace(:,:,elec),.2,'g',linspace(begin_idx, end_idx,size(cond_3_mean_trace,2) ),[])
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
%legend({'inseq +', 'SEM', 'otseq -', 'SEM', 'otseq +', 'SEM', 'inseq -', 'SEM'},'FontSize',10, 'FontWeight', 'bold')
ylim ([mn mx])

%% plot spatial gamma traces
mx = .6
mn = -.47
begin_idx = -.5
end_idx = .5

figure
subplot(1,4,1)
hold on
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_1_mean_trace(:,:,1),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_1_mean_trace(:,:,2),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_1_mean_trace(:,:,3),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_1_mean_trace(:,:,4),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'chan 1' 'chan 2' 'chan 3' 'chan 4'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
title('Inseq +', 'FontSize', 14, 'FontWeight', 'bold')
ylim ([mn mx])

subplot(1,4,3)
hold on
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_2_mean_trace(:,:,1),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_2_mean_trace(:,:,2),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_2_mean_trace(:,:,3),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_2_mean_trace(:,:,4),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'chan 1' 'chan 2' 'chan 3' 'chan 4'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
title('Outseq +', 'FontSize', 14, 'FontWeight', 'bold')
ylim ([mn mx])

subplot(1,4,4)
hold on
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_3_mean_trace(:,:,1),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_3_mean_trace(:,:,2),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_3_mean_trace(:,:,3),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_3_mean_trace(:,:,4),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'chan 1' 'chan 2' 'chan 3' 'chan 4'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
title('InSeq -', 'FontSize', 14, 'FontWeight', 'bold')
ylim ([mn mx])


subplot(1,4,2)
hold on
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_4_mean_trace(:,:,1),1),ones(1,win)/win, 'same'), 'r', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_4_mean_trace(:,:,2),1),ones(1,win)/win, 'same'), 'm', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_4_mean_trace(:,:,3),1),ones(1,win)/win, 'same'), 'b', 'LineWidth', 1.8)
plot(linspace(begin_idx, end_idx,size(cond_1_mean_trace,2) ), conv(nanmean(cond_4_mean_trace(:,:,4),1),ones(1,win)/win, 'same'), 'g', 'LineWidth', 1.8)
legend({'chan 1' 'chan 2' 'chan 3' 'chan 4'}, 'FontSize', 15, 'FontWeight', 'bold')
plot(zeros(1,length(mn:.05:mx)), mn:.05:mx,'Color','k','LineWidth',1)
title('Outseq -', 'FontSize', 14, 'FontWeight', 'bold')
ylim ([mn mx])

