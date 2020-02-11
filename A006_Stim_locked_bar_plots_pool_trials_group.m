clear all;close all;clc
lock = 'p_i' %'p_o' % p_i: poke-in % p_o: poke-out
task = 'welltrained' %novel1, novel2, welltrained
pool_trials = ''
anim_list   = [1:5];
anim_length = length(anim_list);
if length(anim_list)==1
    plot_type = ['indiv_anim' num2str(anim_list)];
else
    plot_type = 'group';
end
cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')

% wavelet params
fs = 1000;
dt = 1/fs;
NumVoices = 32;
a0 = 2^(1/NumVoices);
wavCenterFreq = 6/(2*pi);
minfreq = 3;
maxfreq = 250;
minscale = wavCenterFreq/(maxfreq*dt);
maxscale = wavCenterFreq/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;
freq = wavCenterFreq./(fs*scales.*dt);


cond1_ch1_temp = cell(1,anim_length);
cond1_ch2_temp = cell(1,anim_length);
cond1_ch3_temp = cell(1,anim_length);
cond1_ch4_temp = cell(1,anim_length);

cond2_ch1_temp = cell(1,anim_length);
cond2_ch2_temp = cell(1,anim_length);
cond2_ch3_temp = cell(1,anim_length);
cond2_ch4_temp = cell(1,anim_length);

cond3_ch1_temp = cell(1,anim_length);
cond3_ch2_temp = cell(1,anim_length);
cond3_ch3_temp = cell(1,anim_length);
cond3_ch4_temp = cell(1,anim_length);

cond4_ch1_temp = cell(1,anim_length);
cond4_ch2_temp = cell(1,anim_length);
cond4_ch3_temp = cell(1,anim_length);
cond4_ch4_temp = cell(1,anim_length);

% group all the trail together
for anim_cntr = 1:anim_length
    load(['anim' num2str(anim_list(anim_cntr)) '_spectrogram_data_' lock '_' task])
    
    % chan1
    cond1_ch1_temp{anim_cntr} = norm_freq_acrs_chan_cond_1(:,:,:,1);
    cond2_ch1_temp{anim_cntr} = norm_freq_acrs_chan_cond_2(:,:,:,1);
    cond3_ch1_temp{anim_cntr} = norm_freq_acrs_chan_cond_3(:,:,:,1);
    cond4_ch1_temp{anim_cntr} = norm_freq_acrs_chan_cond_4(:,:,:,1);
    
    % chan2
    cond1_ch2_temp{anim_cntr} = norm_freq_acrs_chan_cond_1(:,:,:,2);
    cond2_ch2_temp{anim_cntr} = norm_freq_acrs_chan_cond_2(:,:,:,2);
    cond3_ch2_temp{anim_cntr} = norm_freq_acrs_chan_cond_3(:,:,:,2);
    cond4_ch2_temp{anim_cntr} = norm_freq_acrs_chan_cond_4(:,:,:,2);
    
    % chan3
    cond1_ch3_temp{anim_cntr} = norm_freq_acrs_chan_cond_1(:,:,:,3);
    cond2_ch3_temp{anim_cntr} = norm_freq_acrs_chan_cond_2(:,:,:,3);
    cond3_ch3_temp{anim_cntr} = norm_freq_acrs_chan_cond_3(:,:,:,3);
    cond4_ch3_temp{anim_cntr} = norm_freq_acrs_chan_cond_4(:,:,:,3);
    
    % chan4
    cond1_ch4_temp{anim_cntr} = norm_freq_acrs_chan_cond_1(:,:,:,4);
    cond2_ch4_temp{anim_cntr} = norm_freq_acrs_chan_cond_2(:,:,:,4);
    cond3_ch4_temp{anim_cntr} = norm_freq_acrs_chan_cond_3(:,:,:,4);
    cond4_ch4_temp{anim_cntr} = norm_freq_acrs_chan_cond_4(:,:,:,4);
end

if strcmp('yes',pool_trials)
    
    % chan1 all conds
    cond1_ch1_temp_temp = cat(3,cond1_ch1_temp{:});
    cond2_ch1_temp_temp = cat(3,cond2_ch1_temp{:});
    cond3_ch1_temp_temp = cat(3,cond3_ch1_temp{:});
    cond4_ch1_temp_temp = cat(3,cond4_ch1_temp{:});
    
    % chan2 all conds
    cond1_ch2_temp_temp = cat(3,cond1_ch2_temp{:});
    cond2_ch2_temp_temp = cat(3,cond2_ch2_temp{:});
    cond3_ch2_temp_temp = cat(3,cond3_ch2_temp{:});
    cond4_ch2_temp_temp = cat(3,cond4_ch2_temp{:});
    
    % chan3 all conds
    cond1_ch3_temp_temp = cat(3,cond1_ch3_temp{:});
    cond2_ch3_temp_temp = cat(3,cond2_ch3_temp{:});
    cond3_ch3_temp_temp = cat(3,cond3_ch3_temp{:});
    cond4_ch3_temp_temp = cat(3,cond4_ch3_temp{:});
    
    % chan4 all conds
    cond1_ch4_temp_temp = cat(3,cond1_ch4_temp{:});
    cond2_ch4_temp_temp = cat(3,cond2_ch4_temp{:});
    cond3_ch4_temp_temp = cat(3,cond3_ch4_temp{:});
    cond4_ch4_temp_temp = cat(3,cond4_ch4_temp{:});
    
    
else
    for anim = 1:anim_length
        % chan1 all conds
        cond1_ch1_temp_temp(:,:,anim) = nanmean(cond1_ch1_temp{anim},3);
        cond2_ch1_temp_temp(:,:,anim) = nanmean(cond2_ch1_temp{anim},3);
        cond3_ch1_temp_temp(:,:,anim) = nanmean(cond3_ch1_temp{anim},3);
        cond4_ch1_temp_temp(:,:,anim) = nanmean(cond4_ch1_temp{anim},3);
        
        % chan2 all conds
        cond1_ch2_temp_temp(:,:,anim) = nanmean(cond1_ch2_temp{anim},3);
        cond2_ch2_temp_temp(:,:,anim) = nanmean(cond2_ch2_temp{anim},3);
        cond3_ch2_temp_temp(:,:,anim) = nanmean(cond3_ch2_temp{anim},3);
        cond4_ch2_temp_temp(:,:,anim) = nanmean(cond4_ch2_temp{anim},3);
        
        % chan3 all conds
        cond1_ch3_temp_temp(:,:,anim) = nanmean(cond1_ch3_temp{anim},3);
        cond2_ch3_temp_temp(:,:,anim) = nanmean(cond2_ch3_temp{anim},3);
        cond3_ch3_temp_temp(:,:,anim)= nanmean(cond3_ch3_temp{anim},3);
        cond4_ch3_temp_temp(:,:,anim) = nanmean(cond4_ch3_temp{anim},3);
        
        % chan4 all conds
        cond1_ch4_temp_temp(:,:,anim) = nanmean(cond1_ch4_temp{anim},3);
        cond2_ch4_temp_temp(:,:,anim) = nanmean(cond2_ch4_temp{anim},3);
        cond3_ch4_temp_temp(:,:,anim) = nanmean(cond3_ch4_temp{anim},3);
        cond4_ch4_temp_temp(:,:,anim) = nanmean(cond4_ch4_temp{anim},3);
    end
end


% now put each of the 4 direc chans in their respective cond matrices
cond_1_all_animals(:,:,:,1) =  cond1_ch1_temp_temp; % freqXtimeXtrialsXchan cond 1
cond_1_all_animals(:,:,:,2) =  cond1_ch2_temp_temp; % freqXtimeXtrialsXchan cond 1
cond_1_all_animals(:,:,:,3) =  cond1_ch3_temp_temp; % freqXtimeXtrialsXchan cond 1
cond_1_all_animals(:,:,:,4) =  cond1_ch4_temp_temp; % freqXtimeXtrialsXchan cond 1

cond_2_all_animals(:,:,:,1) =  cond2_ch1_temp_temp;
cond_2_all_animals(:,:,:,2) =  cond2_ch2_temp_temp;
cond_2_all_animals(:,:,:,3) =  cond2_ch3_temp_temp;
cond_2_all_animals(:,:,:,4) =  cond2_ch4_temp_temp;

cond_3_all_animals(:,:,:,1) =  cond3_ch1_temp_temp;
cond_3_all_animals(:,:,:,2) =  cond3_ch2_temp_temp;
cond_3_all_animals(:,:,:,3) =  cond3_ch3_temp_temp;
cond_3_all_animals(:,:,:,4) =  cond3_ch4_temp_temp;

cond_4_all_animals(:,:,:,1) =  cond4_ch1_temp_temp;
cond_4_all_animals(:,:,:,2) =  cond4_ch2_temp_temp;
cond_4_all_animals(:,:,:,3) =  cond4_ch3_temp_temp;
cond_4_all_animals(:,:,:,4) =  cond4_ch4_temp_temp;

% group spectrograms across axis
edge      = 100;
tickmarks = 1:20:length(freq);
mx        = .75;
mn        = -.75;

if strcmp('p_i', lock)
    pre_stim  = -.5;
    post_stim = 1.5;
    
elseif strcmp('p_o', lock)
    pre_stim  = -.5;
    post_stim = .5;
end

figure;
i = 0;
for i = 1:4
    subplot(1,4,i)
    imagesc(linspace(pre_stim+edge/fs,post_stim-edge/fs,size(cond1_ch1_temp_temp,2)), 1:length(freq), nanmean(cond_1_all_animals(:,edge+1:end-edge,:,i),3))
    colorbar
    colormap jet
    caxis([mn mx])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),  'FontSize', 16, 'FontWeight', 'bold')
    clear mtx
end


%% (fig 2) Inseq+ along axis, spectrograms, one elec at a time
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig2_odorVsrun
%hgexport(gcf, ['group_AL_PM_' task '_' lock '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
i = 0;
for i = 1
    figure
    imagesc(linspace(pre_stim+edge/fs,post_stim-edge/fs,size(cond1_ch1_temp_temp,2)), 1:length(freq), nanmean(cond_1_all_animals(:,edge+1:end-edge,:,i),3))
    colorbar
    colormap jet
    caxis([mn mx])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    %set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),  'FontSize', 16, 'FontWeight', 'bold')
    set(gca,'xtick',[],'ytick',[])

    clear mtx
    hgexport(gcf, [plot_type '_AL_PM_odorInSeqCor_elec_' num2str(i) '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
    hgexport(gcf, [plot_type '_AL_PM_odorInSeqCor_elec_' num2str(i) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
end

%% (fig 2) Inseq+ along axis, barplots and traces: beta
cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')

range = 'beta'
if strcmp('beta',range)
freq_range = freq>20 & freq<40;
colorList    = {'r-o', 'c-o', 'g-o', 'y--o', 'm-o',};
colorListSEM = {'r', 'c', 'g', '--y', 'm'};
elseif strcmp('theta',range)
freq_range = freq>4 & freq<12;
colorList    = {'r--o', 'c-o', 'g-o', 'y-o', 'm-o',};
colorListSEM = {'--r', 'c', 'g', 'y', 'm'};
end
mn = -.8;
mx = .5;
start_time = 610;
end_time   = 1700;

figure;hold on
chan1 = squeeze(nanmean(squeeze(nanmean(cond_1_all_animals(freq_range, start_time:end_time,:,1),1)),1));
chan2 = squeeze(nanmean(squeeze(nanmean(cond_1_all_animals(freq_range, start_time:end_time,:,2),1)),1));
chan3 = squeeze(nanmean(squeeze(nanmean(cond_1_all_animals(freq_range, start_time:end_time,:,3),1)),1));
chan4 = squeeze(nanmean(squeeze(nanmean(cond_1_all_animals(freq_range, start_time:end_time,:,4),1)),1));
bar_vectora1    = [nanmean(chan1) nanmean(chan2) nanmean(chan3) nanmean(chan4)];
bar_vectora1sem = [nanstd(chan1)/sqrt(length(chan1)) nanstd(chan2)/sqrt(length(chan2))...
    nanstd(chan3)/sqrt(length(chan3)) nanstd(chan4)/sqrt(length(chan4))];
b1 = bar(bar_vectora1,'k')
errorbar(1:4,bar_vectora1,bar_vectora1sem,'k*')

colors       = {'r', 'c', 'g', 'y', 'm'};
reps = 1000;
for anim = 1:5
    load(['anim' num2str(anim) '_spectrogram_data_p_i_welltrained.mat'])
    cond_1_mean = [];
    for elec = 1:4
        cond_1_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,elec),1)),1));
    end 
    % traces for indiv animals
    plot(nanmean(cond_1_mean,2) , colorList{anim},'MarkerFaceColor', colors{anim})
    bar_vectora1 = nanmean(cond_1_mean,2);
    bar_vectora1sem = [nanstd(cond_1_mean,0,2)/sqrt(size(cond_1_mean,2)) ];
    errorbar(1:4,bar_vectora1,bar_vectora1sem,colorListSEM{anim})   
    %p = anova1(cond_1_mean')
  [pvals,p_masked,p_fdr] = ElecsAlongAxisStats( reps,cond_1_mean )
end
ylim([mn mx])
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig2_odorVsrun
hgexport(gcf, ['beta_indiv_anim_traces.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
hgexport(gcf, ['beta_indiv_anim_traces.eps'], hgexport('factorystyle'), 'Format', 'eps')


%% save (fig 3) learning spectrogram
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig3_learning
%hgexport(gcf, ['group_AL_PM_' task '_' lock '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')

i = 0;
for i =4
    figure
    imagesc(linspace(pre_stim+edge/fs,post_stim-edge/fs,size(cond1_ch1_temp_temp,2)), 1:length(freq), nanmean(cond_1_all_animals(:,edge+1:end-edge,:,i),3))
    colorbar
    colormap jet
    caxis([mn mx])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),  'FontSize', 16, 'FontWeight', 'bold')
    set(gca,'xtick',[],'ytick',[])

    clear mtx
    hgexport(gcf, [plot_type '_AL_PM_odorCor_'  task '_' lock '_elec_' num2str(i) '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')
    hgexport(gcf, [plot_type '_AL_PM_odorCor_'  task '_' lock '_elec_' num2str(i) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
end


%% save (fig3) traces and stats for indiv animals-learning - single elec
cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')

freq_range = freq>20 & freq<40;  % experimental cond
start_time = 250;
end_time   = 500;

%init for bar plots
cond1_all_anims = zeros(4,1);
cond2_all_anims = zeros(4,1);
cond3_all_anims = zeros(4,1);

% group mean bar plots
for anim = 1:4
    load(['anim' num2str(anim) '_spectrogram_data_p_o_novel1.mat'])
    novel1 = norm_freq_acrs_chan_cond_1;
    
    load(['anim' num2str(anim) '_spectrogram_data_p_o_novel2.mat'])
    novel2 = norm_freq_acrs_chan_cond_1;
    
    load(['anim' num2str(anim) '_spectrogram_data_p_o_welltrained.mat'])
    welltrained = norm_freq_acrs_chan_cond_1;
    cond_1_mean = [];
    cond_2_mean = [];
    cond_3_mean = [];
    
    for elec = 1:4
        cond_1_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(novel1(freq_range, start_time:end_time,:,elec),1)),1));
        cond_2_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(novel2(freq_range, start_time:end_time,:,elec),1)),1));
        cond_3_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(welltrained(freq_range, start_time:end_time,:,elec),1)),1));
    end
    
    cond1 = nanmean(cond_1_mean,1);
    cond2 = nanmean(cond_2_mean,1);
    cond3 = nanmean(cond_3_mean,1);
    
    cond1_all_anims(anim) = nanmean(cond1);
    cond2_all_anims(anim) = nanmean(cond2);
    cond3_all_anims(anim) = nanmean(cond3);
    
end
figure;hold on
bar_vectora1    = [nanmean(cond1_all_anims) nanmean(cond2_all_anims) nanmean(cond3_all_anims) ];
bar_vectora1sem = [nanstd(cond1_all_anims)/sqrt(length(cond1_all_anims)) nanstd(cond2_all_anims)/sqrt(length(cond2_all_anims))...
    nanstd(cond3_all_anims)/sqrt(length(cond3_all_anims))];
b1 = bar(bar_vectora1,'k')
errorbar(1:3,bar_vectora1,bar_vectora1sem,'k*')

% indiv anim traces
colorList = {'r-o', 'c--o', 'g--o', 'y-o'};
colorListSEM = {'-r', '--c', '--g', '-y'};
colors = {'r', 'c', 'g', 'y'};
band_value = '20-40 hz';
for anim = 1:4
    load(['anim' num2str(anim) '_spectrogram_data_p_o_novel1.mat'])
    novel1 = norm_freq_acrs_chan_cond_1;
    
    load(['anim' num2str(anim) '_spectrogram_data_p_o_novel2.mat'])
    novel2 = norm_freq_acrs_chan_cond_1;
    
    load(['anim' num2str(anim) '_spectrogram_data_p_o_welltrained.mat'])
    welltrained = norm_freq_acrs_chan_cond_1;
    cond_1_mean = [];
    cond_2_mean = [];
    cond_3_mean = [];
    cond_1_mean = [];
    cond_2_mean = [];
    cond_3_mean = [];
    
    for elec = 1:4
        cond_1_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(novel1(freq_range, start_time:end_time,:,elec),1)),1));
        cond_2_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(novel2(freq_range, start_time:end_time,:,elec),1)),1));
        cond_3_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(welltrained(freq_range, start_time:end_time,:,elec),1)),1));
    end
    
    cond1 = nanmean(cond_1_mean,1);
    cond2 = nanmean(cond_2_mean,1);
    cond3 = nanmean(cond_3_mean,1);
    
    % traces for indiv animals
    plot([nanmean(cond1) nanmean(cond2) nanmean(cond3)], colorList{anim},'MarkerFaceColor', colors{anim})
    bar_vectora1    = [nanmean(cond1) nanmean(cond2) nanmean(cond3)];
    bar_vectora1sem = [nanstd(cond1)/sqrt(sum(~isnan(cond1))) nanstd(cond2)/sqrt(sum(~isnan(cond2)))  nanstd(cond3)/sqrt(sum(~isnan(cond3)))];
    errorbar(1:3,bar_vectora1,bar_vectora1sem,colorListSEM{anim},'LineWidth',2)%ylim([mn mx])
    set(gca, 'XTick', 1:3, 'XTickLabel', {'novel1' 'novel2' 'welltrained' },'XTickLabelRotation',45)
    
    % unbalanced one way anova
    %  group= zeros( length([cond1'; cond2'; cond3']), 1);
    %  group(1:length(cond1)) = 1;
    %  group(length(cond1)+1:length(cond1)+length(cond2)) = 2;
    %  group(length(cond1)+length(cond2)+1:end) = 3;
    %  p = anova1([cond1'; cond2'; cond3'], group)

    %[ pvals,p_fdr,p_masked ] = LearningStats( 1000,cond1, cond2, cond3 )
end
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end

cd ./fig3_learning
hgexport(gcf, ['group_learning_single_anim_traces_single_elec.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, ['group_learning_single_anim_traces_single_elec.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')

%% save (fig3) traces for indiv animals- beta along proximodistal axis- axis

figure
hold on
%init for bar plots
cond1_all_anims = zeros(5,1);
cond2_all_anims = zeros(5,1);
cond3_all_anims = zeros(5,1);
cond4_all_anims = zeros(5,1);

for anim =1:5
    cond1 = [];
    cond2 = [];
    cond3 = [];
    cond4 = [];
    
    cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
    load(['anim' num2str(anim) '_spectrogram_data_p_o_welltrained.mat'])
    cond1 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,1),1)),1));
    cond2 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,2),1)),1));
    cond3 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,3),1)),1));
    cond4 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,4),1)),1));
    
    % traces for indiv animals
    cond1_all_anims(anim) = nanmean(cond1);
    cond2_all_anims(anim) = nanmean(cond2);
    cond3_all_anims(anim) = nanmean(cond3);
    cond4_all_anims(anim) = nanmean(cond4);

end
figure;hold on
bar_vectora1    = [nanmean(cond1_all_anims) nanmean(cond2_all_anims) nanmean(cond3_all_anims) nanmean(cond4_all_anims)];
bar_vectora1sem = [nanstd(cond1_all_anims)/sqrt(length(cond1_all_anims)) nanstd(cond2_all_anims)/sqrt(length(cond2_all_anims))...
    nanstd(cond3_all_anims)/sqrt(length(cond3_all_anims))  nanstd(cond4_all_anims)/sqrt(length(cond4_all_anims))];
b1 = bar(bar_vectora1,'k')
errorbar(1:4,bar_vectora1,bar_vectora1sem,'k*')


colorList    = {'r-o', 'c-o', 'g-o', 'y-o', 'm--o',};
colorListSEM = {'r', 'c', 'g', 'y', 'm--'};
colors       = {'r', 'c', 'g', 'y', 'm'};
for anim =1:5
    cond1 = [];
    cond2 = [];
    cond3 = [];
    cond4 = [];
    
    cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
    load(['anim' num2str(anim) '_spectrogram_data_p_o_welltrained.mat'])
    cond1 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,1),1)),1));
    cond2 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,2),1)),1));
    cond3 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,3),1)),1));
    cond4 = squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,4),1)),1));
    
    % traces for indiv animals
    plot([nanmean(cond1) nanmean(cond2) nanmean(cond3) nanmean(cond4)], colorList{anim},'MarkerFaceColor', colors{anim})
    bar_vectora1 = [nanmean(cond1) nanmean(cond2) nanmean(cond3) nanmean(cond4)];
    bar_vectora1sem = [nanstd(cond1)/sqrt(sum(~isnan(cond1))) nanstd(cond2)/sqrt(sum(~isnan(cond2)))...
        nanstd(cond3)/sqrt(sum(~isnan(cond3))) nanstd(cond4)/sqrt(sum(~isnan(cond4)))];
    errorbar(1:4,bar_vectora1,bar_vectora1sem, colorListSEM{anim})%ylim([mn mx])
    set(gca, 'XTick', 1:4, 'XTickLabel', {'proimal CA1' '' '' 'distal CA1' },'XTickLabelRotation',45)
    
    % run anova then perm testing
%   p = anova1([cond1' cond2' cond3' cond4'])
  % [pvals,p_masked,p_fdr] = ElecsAlongAxisStats( 1000,[cond1' cond2' cond3' cond4'] )
end


if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig3_learning
hgexport(gcf, ['group_betapower_axis_single_anim_traces' lock '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, ['group_betapower_axis_single_anim_traces' lock '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')

%% fig 4,sup - main bar plots - all 4 trials types - single elec of average elecs across axis

fig_num = 'supb' % 4, supa, supb, supc

if strcmp('4',fig_num)
    % fig 4 params
    band_value = '20-40 hz';
    freq_range = freq>20 & freq<40;
    if strcmp('p_o', lock)
        start_time = 250;
        end_time   = 500;
        mn = -.3
        mx = .8
    elseif  strcmp('p_i', lock)
        start_time = 610; %110ms
        end_time   = 900; %400ms
    end
    
elseif strcmp('supa',fig_num) % slow gamma
    % fig 6a params
    band_value = '25-55 hz';
    freq_range = freq>25 & freq<55;
    start_time = 250; %
    end_time   = 500; %

elseif strcmp('supb',fig_num) % fast gamma
    % fig 6a params
    band_value = '60-100 hz';
    freq_range = freq>60 & freq<100;
    start_time = 610; %
    end_time   = 810; %

elseif strcmp('supc',fig_num) % theta
    % fig 6a params
    band_value = '4-8 hz';
    freq_range = freq>4 & freq<8;
    start_time = 610; %
    end_time   = 900; %

end
if strcmp('p_i', lock)
    mn = -.8
    mx = .5
end
cond_1_all_animals_plot = nanmean(cond_1_all_animals,4);
cond_2_all_animals_plot = nanmean(cond_2_all_animals,4);
cond_3_all_animals_plot = nanmean(cond_3_all_animals,4);
cond_4_all_animals_plot = nanmean(cond_4_all_animals,4);

cond_1_mean = squeeze(nanmean(squeeze(nanmean(cond_1_all_animals_plot(freq_range, start_time:end_time,:),1)),1)); % pooled trialsXchans
cond_1_ntrls = length(cond_1_mean) - sum(isnan(cond_1_mean));

cond_2_mean = squeeze(nanmean(squeeze(nanmean(cond_2_all_animals_plot(freq_range, start_time:end_time,:),1)),1)); % pooled trialsXchans
cond_2_ntrls = length(cond_2_mean) - sum(isnan(cond_2_mean));

cond_3_mean = squeeze(nanmean(squeeze(nanmean(cond_3_all_animals_plot(freq_range, start_time:end_time,:),1)),1)); % pooled trialsXchans
cond_3_ntrls = length(cond_3_mean) - sum(isnan(cond_3_mean));

cond_4_mean = squeeze(nanmean(squeeze(nanmean(cond_4_all_animals_plot(freq_range, start_time:end_time,:),1)),1)); % pooled trialsXchans
cond_4_ntrls = length(cond_4_mean) - sum(isnan(cond_4_mean));

figure
bar_vectora1 = [nanmean(cond_1_mean) nanmean(cond_4_mean) nanmean(cond_2_mean) nanmean(cond_3_mean)  ];
bar(bar_vectora1, 'k')
hold on
bar_vectora1sem = [nanstd(cond_1_mean)/sqrt(cond_1_ntrls) nanstd(cond_4_mean)/sqrt(cond_4_ntrls) nanstd(cond_2_mean)/sqrt(cond_2_ntrls) nanstd(cond_3_mean)/sqrt(cond_3_ntrls)  ];
errorbar(1:4,bar_vectora1,bar_vectora1sem, 'k')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'outseq -' 'outseq +' 'inseq -' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
%ylim([mn mx])


% draw indiv anim data
colorList    = {'ro', 'co', 'go', 'yo', 'mo',};
colorListSEM = {' r', ' c', ' g', ' y', ' m'};
colors = {'r', 'c', 'g', 'y', 'm'};

hold on
for anim =1:5
    cond1 = [];
    cond2 = [];
    cond3 = [];
    cond4 = [];
    
    cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
    load(['anim' num2str(anim) '_spectrogram_data_' lock '_welltrained.mat'])
    cond1 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,:),1)),1)),2);
    cond2 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_2(freq_range, start_time:end_time,:,:),1)),1)),2);
    cond3 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_3(freq_range, start_time:end_time,:,:),1)),1)),2);
    cond4 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_4(freq_range, start_time:end_time,:,:),1)),1)),2);
    
    % traces for indiv animals
    plot([nanmean(cond1)  nanmean(cond4) nanmean(cond2) nanmean(cond3)], colorList{anim},'MarkerFaceColor', colors{anim})
    bar_vectora1 = [nanmean(cond1) nanmean(cond4) nanmean(cond2) nanmean(cond3)];
    bar_vectora1sem = [nanstd(cond1)/sqrt(sum(~isnan(cond1))) nanstd(cond4)/sqrt(sum(~isnan(cond4)))...
        nanstd(cond2)/sqrt(sum(~isnan(cond2))) nanstd(cond3)/sqrt(sum(~isnan(cond3))) ];
    errorbar(1:4,bar_vectora1,bar_vectora1sem, colorListSEM{anim},'LineStyle','none')
    set(gca, 'XTick', 1:4, 'XTickLabel', {'InSeq+' 'OutSeq-' 'OutSeq+' 'InSeq-' },'XTickLabelRotation',45)

end

if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end

if strcmp('4',fig_num)
    cd ./fig4_beta
    hgexport(gcf, ['group_betapower_4TrialTypes_single_anim_traces_' lock '.eps'], hgexport('factorystyle'), 'Format', 'eps')
    hgexport(gcf, ['group_betapower_4TrialTypes_single_anim_traces_' lock '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
else
    cd ./sup
    hgexport(gcf, ['group_gamma_' band_value '_4TrialTypes_single_anim_traces_' lock '_' num2str(start_time) '_' num2str(end_time) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
    hgexport(gcf, ['group_gamma_' band_value '_4TrialTypes_single_anim_traces_' lock '_' num2str(start_time) '_' num2str(end_time) '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')  
end
%% fig 4,sup - match or corr  plots - single elec of average elecs across axis + indiv anim traces
plot_cond = 'correct'
if strcmp('p_o', lock)
    mx = .6;
    mn = -.1;
    if strcmp('match', plot_cond)
        colorList    = {'--ro', '-co', '--go', '--yo', '-mo',};
        colorListSEM = {' --r', ' -c', '--g', '--y', '-m'};
        
    else strcmp('correct', plot_cond)
        colorList    = {'--ro', '--co', '--go', '--yo', '-mo',};
        colorListSEM = {' --r', ' --c', '--g', '--y', '-m'};
    end
elseif  strcmp('p_i', lock)
    if strcmp('4',fig_num)
        colorList    = {'--ro', '--co', '--go', '--yo', '--mo',};
        colorListSEM = {'--r', '--c', '--g', '--y', '--m'};
    elseif strcmp('supa',fig_num)
        if strcmp('match', plot_cond)
            colorList    = {'--ro', '--co', '--go', '--yo', '--mo',};
            colorListSEM = {' --r', ' --c', '--g', '--y', '--m'};
            
        else strcmp('correct', plot_cond)
            colorList    = {'--ro', '--co', '--go', '--yo', '-mo',};
            colorListSEM = {' --r', ' --c', '--g', '--y', '-m'};
        end
    elseif strcmp('supb',fig_num)
        if strcmp('match', plot_cond)
            colorList    = {'--ro', '--co', '--go', '--yo', '-mo',};
            colorListSEM = {' --r', ' --c', '--g', '--y', '-m'};
        else strcmp('correct', plot_cond)
            colorList    = {'--ro', '--co', '--go', '--yo', '-mo',};
            colorListSEM = {' --r', ' --c', '--g', '--y', '-m'};
        end
        
    elseif strcmp('supc',fig_num)
        if strcmp('match', plot_cond)
            colorList    = {'--ro', '--co', '--go', '--yo', '-mo',};
            colorListSEM = {' --r', ' --c', '--g', '--y', '-m'};
            
        elseif strcmp('correct', plot_cond)
            colorList    = {'--ro', '--co', '--go', '--yo', '--mo'};
            colorListSEM = {' --r', ' --c', '--g', '--y', '--m'};
        end
    end

end
if strcmp('match', plot_cond)
    cond1 = nanmean([cond_1_mean; cond_4_mean],1);
    cond2 = nanmean([cond_2_mean; cond_3_mean],1);
    legend_entry = {'match', 'mismatch'};
else strcmp('correct', plot_cond)
    cond1 = nanmean([cond_1_mean; cond_2_mean],1);
    cond2 = nanmean([cond_3_mean; cond_4_mean],1);
    legend_entry = {'corr', 'incorr'};
end

figure
bar_vectora1 = [nanmean(cond1) nanmean(cond2)];
bar(bar_vectora1, 'k')
hold on
bar_vectora1sem = [nanstd(cond1)/sqrt(sum(~isnan(cond1))) nanstd(cond2)/sqrt(sum(~isnan(cond2)))];
errorbar(1:2, bar_vectora1,bar_vectora1sem, 'k')%ylim([mn mx])
set(gca, 'XTick', 1:2, 'XTickLabel', legend_entry,'XTickLabelRotation',45)
%set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])
colors = {'r', 'c', 'g', 'y', 'm'};

hold on
for anim =1:5
    cond1 = [];
    cond2 = [];
    cond3 = [];
    cond4 = [];
    
    cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
    load(['anim' num2str(anim) '_spectrogram_data_' lock '_welltrained.mat'])
    cond1 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(freq_range, start_time:end_time,:,:),1)),1)),2);
    cond2 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_2(freq_range, start_time:end_time,:,:),1)),1)),2);
    cond3 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_3(freq_range, start_time:end_time,:,:),1)),1)),2);
    cond4 = nanmean(squeeze(nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_4(freq_range, start_time:end_time,:,:),1)),1)),2);
    
    if strcmp('match', plot_cond)
        cond1_plt = [cond1' cond4'];
        cond2_plt = [cond2' cond3'];
    else strcmp('correct', plot_cond)
        cond1_plt = [cond1' cond2'];
        cond2_plt =  [cond3' cond4'];
    end
    % traces for indiv animals
    plot([nanmean(cond1_plt) nanmean(cond2_plt)], colorList{anim},'MarkerFaceColor', colors{anim})
    bar_vectora1 = [nanmean(cond1_plt) nanmean(cond2_plt)];
    bar_vectora1sem = [nanstd(cond1_plt)/sqrt(sum(~isnan(cond1_plt))) nanstd(cond2_plt)/sqrt(sum(~isnan(cond2_plt)))]
    errorbar(1:2,bar_vectora1,bar_vectora1sem, colorListSEM{anim})%ylim([mn mx])
    set(gca, 'XTick', 1:2, 'XTickLabel', legend_entry,'XTickLabelRotation',45)

   % run stats per anim
%     reps  = 1000;
%     adata = cond1_plt(~isnan(cond1_plt));
%     bdata = cond2_plt(~isnan(cond2_plt));
%     p1(anim)= permutation_unpaired(adata, bdata, reps)
end


% save
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
if strcmp('4',fig_num)
    cd ./fig4_beta
elseif strcmp('supa',fig_num) || strcmp('supb',fig_num)
    cd ./sup
end
hgexport(gcf, ['group_' band_value '_' plot_cond '_single_anim_traces_' lock '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, ['group_' band_value '_' plot_cond '_single_anim_traces_' lock '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')

%% fig 5 - cluster - single value per anim
cd('D:\Gattas\ephys_data_final\group_plots\cluster_matrices')
zmap_index = zeros(size(norm_freq_acrs_chan_cond_1(:,:,1)));
load('low_freq_group_corrected_zmap_p_i.mat')
zmap_index(:, 501:1000) = corrected_zmap;

% divide clusters in 3 subclusters
full_cluster = zmap_index;

zmap_theta_subset1      = zmap_index;
zmap_theta_subset1(freq>9,:) = 0;

zmap_subset2 = zmap_index;
zmap_subset2(freq<9,:)  = 0;
zmap_subset2(freq>12.5,:) = 0;

zmap_theta_harm_subset3= zmap_index;
zmap_theta_harm_subset3(freq<12.5,:)=0;

% select cluster range
subset = 'one'
if strcmp('one', subset)
    zmap_index = zmap_theta_subset1;
    mx = 0.8
    mn = -0.7;
elseif strcmp('two', subset)
    zmap_index  = zmap_subset2;
        mx = 0.9
    mn = -1.8;
elseif strcmp('three', subset)
    zmap_index  = zmap_theta_harm_subset3;
        mx = 0.5
    mn = -1;
else
    zmap_index  = full_cluster;
    mx = 0.5;
    mn = -0.5;
end

% first plot spectrogram

figure
imagesc(linspace(pre_stim+edge/fs,post_stim-edge/fs,size(cond1_ch1_temp_temp,2)), 1:length(freq), nanmean(cond_1_all_animals(:,edge+1:end-edge,:,1),3))
hold on;plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
contour(linspace(pre_stim+(edge/fs),post_stim-(edge/fs),size(cond1_ch1_temp_temp,2)), 1:length(freq), zmap_index,40, 'LineColor', 'k', 'LineWidth', 2.5)
set(gca,'YTick',tickmarks, 'YTickLabel',round(freq(tickmarks)));
colorbar;colormap jet;caxis([-0.8 0.8]);xlabel('time');ylabel('freq');title('InSeq+')
% save
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig5_cluster
hgexport(gcf, ['group_InSeqCorr_cluster_subset' subset '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
hgexport(gcf, ['group_InSeqCorr_cluster_subset'  subset '.eps'], hgexport('factorystyle'), 'Format', 'eps')


% bar plots
for anim = 1:5
temp = cond_1_all_animals_plot(:,:,anim); 
cond_1_mean(anim) = mean(mean(temp(logical(zmap_index)),1),1);

temp = []; temp = cond_2_all_animals_plot(:,:,anim); 
cond_2_mean(anim) = mean(mean(temp(logical(zmap_index)),1),1)

temp = []; temp = cond_3_all_animals_plot(:,:,anim); 
cond_3_mean(anim) = mean(mean(temp(logical(zmap_index)),1),1);

temp = []; temp = cond_4_all_animals_plot(:,:,anim); 
cond_4_mean(anim) = mean(mean(temp(logical(zmap_index)),1),1);
end



figure
bar_vectora1 = [nanmean(cond_1_mean) nanmean(cond_4_mean) nanmean(cond_2_mean) nanmean(cond_3_mean)  ];
bar(bar_vectora1, 'b')
hold on
bar_vectora1sem = [nanstd(cond_1_mean)/sqrt(length(cond_1_mean)) nanstd(cond_4_mean)/sqrt(length(cond_4_mean))...
    nanstd(cond_2_mean)/sqrt(length(cond_2_mean)) nanstd(cond_3_mean)/sqrt(length(cond_2_mean))  ];
errorbar(1:4,bar_vectora1,bar_vectora1sem, 'kx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'outseq -' 'outseq +' 'inseq -' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])

% draw indiv anim data
colorList    = {'ro', 'co', 'go', 'yo', 'mo',};
colorListSEM = {' r', ' c', ' g', ' y', ' m'};
colors = {'r', 'c', 'g', 'y', 'm'};

hold on
for anim =1:5
    
    cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
    load(['anim' num2str(anim) '_spectrogram_data_' lock '_welltrained.mat'])
    cond1 = zeros(size(norm_freq_acrs_chan_cond_1,3), size(norm_freq_acrs_chan_cond_1,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_1,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_1,3)
            temp = norm_freq_acrs_chan_cond_1(:,:,iTrl,iElec);
            cond1 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond1 = nanmean(cond1,2);
    
    cond2 = zeros(size(norm_freq_acrs_chan_cond_2,3), size(norm_freq_acrs_chan_cond_2,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_2,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_2,3)
            temp = norm_freq_acrs_chan_cond_2(:,:,iTrl,iElec);
            cond2 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond2 = nanmean(cond2,2);
    
    cond3 = zeros(size(norm_freq_acrs_chan_cond_3,3), size(norm_freq_acrs_chan_cond_3,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_3,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_3,3)
            temp = norm_freq_acrs_chan_cond_3(:,:,iTrl,iElec);
            cond3 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond3 = nanmean(cond3,2);
    
    cond4 = zeros(size(norm_freq_acrs_chan_cond_4,3), size(norm_freq_acrs_chan_cond_4,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_4,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_4,3)
            temp = norm_freq_acrs_chan_cond_4(:,:,iTrl,iElec);
            cond4 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond4 = nanmean(cond4,2);
    % traces for indiv animals
    plot([nanmean(cond1)  nanmean(cond4) nanmean(cond2) nanmean(cond3)], colorList{anim},'MarkerFaceColor', colors{anim},'LineWidth',2,'MarkerSize',10)
    bar_vectora1 = [nanmean(cond1) nanmean(cond4) nanmean(cond2) nanmean(cond3)];
    bar_vectora1sem = [nanstd(cond1)/sqrt(sum(~isnan(cond1))) nanstd(cond4)/sqrt(sum(~isnan(cond4)))...
        nanstd(cond2)/sqrt(sum(~isnan(cond2))) nanstd(cond3)/sqrt(sum(~isnan(cond3))) ];
    errorbar(1:4,bar_vectora1,bar_vectora1sem, colorListSEM{anim},'LineStyle','none')%ylim([mn mx])
    set(gca, 'XTick', 1:4, 'XTickLabel', {'InSeq+' 'OutSeq-' 'OutSeq+' 'InSeq-' },'XTickLabelRotation',45)
end
% save
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig5_cluster
hgexport(gcf, ['group_4TrialTypes_cluster_subset'  subset '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, ['group_4TrialTypes_cluster_subset'  subset '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')

%% match, correct, cluster plots
plot_cond = 'correct'
if strcmp('match', plot_cond)
    cond1 = [cond_1_mean cond_4_mean];
    cond2 = [cond_2_mean cond_3_mean];
    legend_entry = {'match', 'mismatch'};
    
else strcmp('correct', plot_cond)
    cond1 = [cond_1_mean cond_2_mean];
    cond2 = [cond_3_mean cond_4_mean];
    legend_entry = {'corr', 'incorr'};
end

figure
bar_vectora1 = [nanmean(cond1) nanmean(cond2)];
bar(bar_vectora1, 'b')
hold on
bar_vectora1sem = [nanstd(cond1)/sqrt(sum(~isnan(cond1))) nanstd(cond2)/sqrt(sum(~isnan(cond2)))];
errorbar(1:2, bar_vectora1,bar_vectora1sem, 'kx')%ylim([mn mx])
set(gca, 'XTick', 1:2, 'XTickLabel', legend_entry,'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
ylim([mn mx])
colors = {'r', 'c', 'g', 'y', 'm'};


% draw indiv animal traces
hold on
if strcmp('', subset)
    mx = .4;
    mn = -.4;
    if strcmp('match', plot_cond)
        colorList    = {'--ro', '--co', '-go', '--yo', '-mo',};
        colorListSEM = {' --r', ' --c', '-g', '--y', '-m'};
        
    else strcmp('correct', plot_cond)
        colorList    = {'--ro', '--co', '--go', '--yo', '-mo',};
        colorListSEM = {' --r', ' --c', '--g', '--y', '-m'};
    end
elseif strcmp('one', subset)
    mx = .4;
    mn = -.4;
    if strcmp('match', plot_cond)
        colorList    = {'--ro', '--co', '-go', '--yo', '--mo',};
        colorListSEM = {' --r', ' --c', '-g', '--y', '--m'};
        
    else strcmp('correct', plot_cond)
        colorList    = {'--ro', '--co', '--go', '--yo', '-mo',};
        colorListSEM = {' --r', ' --c', '--g', '--y', '-m'};
    end
    
elseif strcmp('two', subset)
    mx = .2;
    mn = -.9;
    colorList    = {'--ro', '--co', '--go', '--yo', '--mo',};
    colorListSEM = {'--r', ' --c', '--g', '--y', '--m'};
elseif strcmp('three', subset)
    if strcmp('match', plot_cond)
        mn = -.4; mx=.2;
            colorList    = {'--ro', '--co', '--go', '--yo', '--mo',};
    colorListSEM = {'--r', ' --c', '--g', '--y', '--m'};
    else strcmp('correct', plot_cond)
        mn = -.4; mx=.4;
            colorList    = {'--ro', '-co', '--go', '--yo', '-mo',};
    colorListSEM = {'--r', ' -c', '--g', '--y', '-m'};
    end
end

for anim = 1:5
    
    cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
    load(['anim' num2str(anim) '_spectrogram_data_' lock '_welltrained.mat'])
    cond1 = zeros(size(norm_freq_acrs_chan_cond_1,3), size(norm_freq_acrs_chan_cond_1,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_1,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_1,3)
            temp = norm_freq_acrs_chan_cond_1(:,:,iTrl,iElec);
            cond1 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond1 = nanmean(cond1,2);
    
    cond2 = zeros(size(norm_freq_acrs_chan_cond_2,3), size(norm_freq_acrs_chan_cond_2,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_2,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_2,3)
            temp = norm_freq_acrs_chan_cond_2(:,:,iTrl,iElec);
            cond2 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond2 = nanmean(cond2,2);
    
    cond3 = zeros(size(norm_freq_acrs_chan_cond_3,3), size(norm_freq_acrs_chan_cond_3,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_3,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_3,3)
            temp = norm_freq_acrs_chan_cond_3(:,:,iTrl,iElec);
            cond3 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond3 = nanmean(cond3,2);
    
    cond4 = zeros(size(norm_freq_acrs_chan_cond_4,3), size(norm_freq_acrs_chan_cond_4,4));
    for iElec = 1:size(norm_freq_acrs_chan_cond_4,4)
        for iTrl = 1:size(norm_freq_acrs_chan_cond_4,3)
            temp = norm_freq_acrs_chan_cond_4(:,:,iTrl,iElec);
            cond4 (iTrl,iElec) =  nanmean(temp(logical(zmap_index)));
        end
    end
    cond4 = nanmean(cond4,2);
    
    
    if strcmp('match', plot_cond)
        cond1_plt = [cond1' cond4'];
        cond2_plt = [cond2' cond3'];
    else strcmp('correct', plot_cond)
        cond1_plt = [cond1' cond2'];
        cond2_plt =  [cond3' cond4'];
    end
    % traces for indiv animals
    plot([nanmean(cond1_plt) nanmean(cond2_plt)], colorList{anim},'MarkerFaceColor', colors{anim},'LineWidth',2,'MarkerSize',10)
    bar_vectora1 = [nanmean(cond1_plt) nanmean(cond2_plt)];
    bar_vectora1sem = [nanstd(cond1_plt)/sqrt(sum(~isnan(cond1_plt))) nanstd(cond2_plt)/sqrt(sum(~isnan(cond2_plt)))]
    errorbar(1:2,bar_vectora1,bar_vectora1sem, colorListSEM{anim})%ylim([mn mx])
    set(gca, 'XTick', 1:2, 'XTickLabel', legend_entry,'XTickLabelRotation',45)
    ylim([mn mx])
    
    % run stats per anim
    reps  = 1000;
    adata = cond1_plt(~isnan(cond1_plt));
    bdata = cond2_plt(~isnan(cond2_plt));
    p1(anim)= permutation_unpaired(adata, bdata, reps)
end
% save
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig5_cluster
hgexport(gcf, ['group_' plot_cond '_cluster_subset'  subset '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, ['group_' plot_cond '_cluster_subset'  subset '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bar Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig3: Learning bar plots: welltrained, novel 1, novel 2 both across axis and single elec
% load all learning data
cd('D:\Gattas\ephys_data_final\group_plots\spectrogram_matrices_all_conds_per_animal')
if strcmp('yes',pool_trials)
    % load matrices
    load('novel1_InSeq_learning_pooled_trials.mat')% cond1
    load('novel2_InSeq_learning_pooled_trials.mat')% cond2
    load('welltrained_InSeq_learning_pooled_trials.mat')% cond3
    % then cd for fig
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim')
else %single val per anim
    load('novel1_InSeq_learning.mat')% cond1
    load('novel2_InSeq_learning.mat')% cond2
    load('welltrained_InSeq_learning.mat')% cond3
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim')
end
cd ./fig3_learning

band_value = '20-40 hz';
freq_range = freq>20 & freq<40;  % experimental cond
start_time = 250;
end_time   = 500;
figure
mn = 0
mx = .5

cond1_single_elec = [];
cond2_single_elec = [];
cond3_single_elec = [];
for elec = 1:4
    cond_1_mean = squeeze(nanmean(squeeze(nanmean(cond1(freq_range, start_time:end_time,:,elec),1)),1)); % pooled trialsXchans
    cond_1_ntrls = length(cond_1_mean) - sum(isnan(cond_1_mean));
    
    cond_2_mean = squeeze(nanmean(squeeze(nanmean(cond2(freq_range, start_time:end_time,:,elec),1)),1)); % pooled trialsXchans
    cond_2_ntrls = length(cond_2_mean) - sum(isnan(cond_2_mean));
    
    cond_3_mean = squeeze(nanmean(squeeze(nanmean(cond3(freq_range, start_time:end_time,:,elec),1)),1)); % pooled trialsXchans
    cond_3_ntrls = length(cond_3_mean) - sum(isnan(cond_3_mean));
    
    
    subplot(1, 4, elec)
    bar_vectora1 = [nanmean(cond_3_mean) nanmean(cond_2_mean) nanmean(cond_1_mean) ];
    bar(bar_vectora1, 'b')
    hold on
    bar_vectora1sem = [nanstd(cond_3_mean)/sqrt(cond_3_ntrls) nanstd(cond_2_mean)/sqrt(cond_2_ntrls)  nanstd(cond_1_mean)/sqrt(cond_1_ntrls)];
    errorbar(1:3,bar_vectora1,bar_vectora1sem, 'rx')%ylim([mn mx])
    set(gca, 'XTick', 1:3, 'XTickLabel', {'welltrained' 'novel2' 'novel1' },'XTickLabelRotation',45)
    set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    clear bar_vectora1 bar_vectora1std
    ylim([mn mx])
    
    
    hgexport(gcf, ['group_learning_barplots_single_elec.eps'], hgexport('factorystyle'), 'Format', 'eps')
    % single elec
    cond1_single_elec = [cond1_single_elec; cond_1_mean];
    cond2_single_elec = [cond2_single_elec; cond_2_mean];
    cond3_single_elec = [cond3_single_elec; cond_3_mean];
    
end

hgexport(gcf, ['group_learning_barplots_axis.eps'], hgexport('factorystyle'), 'Format', 'eps')

% run stats
% reps  = 1000;
% adata = cond_1_mean;
% bdata = cond_2_mean;
% p1     = permutation_paired(adata, bdata, reps)
%
% adata = cond_1_mean;
% bdata = cond_3_mean;
% p2    = permutation_paired(adata, bdata, reps)
%
% adata = cond_2_mean;
% bdata = cond_3_mean;
% p3    = permutation_paired(adata, bdata, reps)
%
% alpha = .05
% pvals = [p1 p2 p3];
% [p_fdr, p_masked] = fdr( pvals, alpha)
% disp(['novel 2 vs novel1 : ' num2str(pvals(1))])
% disp(['welltrained vs novel1 : ' num2str(pvals(2))])
% disp(['welltrained vs novel2 : ' num2str(pvals(3))])
% p_fdr


% single elec plot
cond1_single_elec = nanmean(cond1_single_elec,1);
cond2_single_elec = nanmean(cond2_single_elec,1);
cond3_single_elec = nanmean(cond3_single_elec,1);
bar_vectora1 = [nanmean(cond3_single_elec) nanmean(cond2_single_elec) nanmean(cond1_single_elec) ];
figure
bar(bar_vectora1, 'b')
hold on
bar_vectora1sem = [nanstd(cond3_single_elec)/sqrt(cond_3_ntrls) nanstd(cond2_single_elec)/sqrt(cond_2_ntrls)  nanstd(cond1_single_elec)/sqrt(cond_1_ntrls)];
errorbar(1:3,bar_vectora1,bar_vectora1sem, 'rx')%ylim([mn mx])
set(gca, 'XTick', 1:3, 'XTickLabel', {'welltrained' 'novel2' 'novel1' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])
hgexport(gcf, ['group_learning_barplots_single_elec.eps'], hgexport('factorystyle'), 'Format', 'eps')

% % run stats
% reps  = 1000;
% adata = cond1_single_elec;
% bdata = cond2_single_elec;
% p1     = permutation_paired(adata, bdata, reps)
%
% adata = cond1_single_elec;
% bdata = cond3_single_elec;
% p2    = permutation_paired(adata, bdata, reps)
%
% adata = cond2_single_elec;
% bdata = cond3_single_elec;
% p3    = permutation_paired(adata, bdata, reps)
%
% alpha = .05
% pvals = [p1 p2 p3];
% [p_fdr, p_masked] = fdr( pvals, alpha)
% disp(['novel 2 vs novel1 : ' num2str(pvals(1))])
% disp(['welltrained vs novel1 : ' num2str(pvals(2))])
% disp(['welltrained vs novel2 : ' num2str(pvals(3))])
% p_fdr




%% Fig 4 and fig 6 (group-pool trials + single val per animal): plot bar plots along axis, barplots for a single elec for a freq range beta, corr incorr, match mistmach
% fig 4
cd ..
cd ./fig4_beta

% % fig 6
% cd ./fig6_gamma

%  plot bar plots along axis
mn = -.28
mx = 0.6
%band_value = '25-55 hz';
%band_value = '60-100 hz';
band_value = '20-40 hz'
start_time = 250;
end_time   = 500;

% start_time = 500;
% end_time   = 750;

if strcmp('25-55 hz',band_value)
    freq_range = freq>25 & freq<55;
elseif strcmp('60-100 hz',band_value)
    freq_range = freq>60 & freq<100;
elseif strcmp('20-40 hz',band_value)
    freq_range = freq>20 & freq<40;
end

% poke out: 250:500
% PI: 610:900 (110-400ms)
% PI: 250:500(-250-0ms)
% PI: 500:750(0-250)
figure
cond1_all_elecs = [];
cond2_all_elecs = [];
cond3_all_elecs = [];
cond4_all_elecs = [];

for elec = 1:4
    cond_1_mean = squeeze(nanmean(squeeze(nanmean(cond_1_all_animals(freq_range, start_time:end_time,:,elec),1)),1)); % pooled trialsXchans
    cond_1_ntrls = length(cond_1_mean) - sum(isnan(cond_1_mean));
    
    cond_2_mean = squeeze(nanmean(squeeze(nanmean(cond_2_all_animals(freq_range, start_time:end_time,:,elec),1)),1)); % pooled trialsXchans
    cond_2_ntrls = length(cond_2_mean) - sum(isnan(cond_2_mean));
    
    cond_3_mean = squeeze(nanmean(squeeze(nanmean(cond_3_all_animals(freq_range, start_time:end_time,:,elec),1)),1)); % pooled trialsXchans
    cond_3_ntrls = length(cond_3_mean) - sum(isnan(cond_3_mean));
    
    cond_4_mean = squeeze(nanmean(squeeze(nanmean(cond_4_all_animals(freq_range, start_time:end_time,:,elec),1)),1)); % pooled trialsXchans
    cond_4_ntrls = length(cond_4_mean) - sum(isnan(cond_4_mean));
    
    subplot(1, 4, elec)
    bar_vectora1 = [nanmean(cond_1_mean) nanmean(cond_4_mean) nanmean(cond_2_mean) nanmean(cond_3_mean)  ];
    bar(bar_vectora1, 'b')
    hold on
    bar_vectora1sem = [nanstd(cond_1_mean)/sqrt(cond_1_ntrls) nanstd(cond_4_mean)/sqrt(cond_4_ntrls) nanstd(cond_2_mean)/sqrt(cond_2_ntrls) nanstd(cond_3_mean)/sqrt(cond_3_ntrls)  ];
    errorbar(1:4,bar_vectora1,bar_vectora1sem, 'kx')%ylim([mn mx])
    set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'outseq -' 'outseq +' 'inseq -' },'XTickLabelRotation',45)
    %set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    clear bar_vectora1 bar_vectora1std
    ylim([mn mx])
    
    cond1_all_elecs = [cond1_all_elecs cond_1_mean];
    cond2_all_elecs = [cond2_all_elecs cond_2_mean];
    cond3_all_elecs = [cond3_all_elecs cond_3_mean];
    cond4_all_elecs = [cond4_all_elecs cond_4_mean];
    
end
suptitle (band_value)

%hgexport(gcf, [ plot_type '_Barplots_Axis_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
%hgexport(gcf, [plot_type '_Barplots_Axis_' band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')
% run stats

% reps  = 1000;
% adata = cond_1_mean;
% bdata = cond_4_mean;
% p1     = permutation_unpaired(adata, bdata, reps)
%
% adata = cond_1_mean;
% bdata = cond_2_mean;
% p2    = permutation_unpaired(adata, bdata, reps)
%
% adata = cond_1_mean;
% bdata = cond_3_mean;
% p3    = permutation_unpaired(adata, bdata, reps)
%
% reps  = 1000;
% adata = cond_2_mean;
% bdata = cond_3_mean;
% p4     = permutation_paired(adata, bdata, reps)
%
% adata = cond_2_mean;
% bdata = cond_3_mean;
% p5    = permutation_unpaired(adata, bdata, reps)
%
% adata = cond_3_mean;
% bdata = cond_4_mean;
% p6    = permutation_unpaired(adata, bdata, reps)
%
%
% alpha = .05
% pvals = [p1 p2 p3 p4 p5 p6];
% [p_fdr, p_masked] = fdr( pvals, alpha)
% disp(['inseq+ vs outseq- : ' num2str(pvals(3))])
% disp(['inseq+ vs outseq+ : ' num2str(pvals(1))])
% disp(['inseq+ vs inseq- : ' num2str(pvals(2))])
% disp(['outseq+ vs outseq- : ' num2str(pvals(5))])
% disp(['inseq- vs outseq- : ' num2str(pvals(6))])
% disp(['outseq+ vs inseq- : ' num2str(pvals(4))])
% p_fdr


% Fig4:one bar plot for entire axis
% get trl num
cond1_all_elecs_ntrls = length(cond1_all_elecs) - sum(isnan(cond1_all_elecs));
cond2_all_elecs_ntrls = length(cond2_all_elecs) - sum(isnan(cond2_all_elecs));
cond3_all_elecs_ntrls = length(cond3_all_elecs) - sum(isnan(cond3_all_elecs));
cond4_all_elecs_ntrls = length(cond4_all_elecs) - sum(isnan(cond4_all_elecs));

figure
mx = .35
mn =0
bar_vectora1 = [nanmean(cond1_all_elecs) nanmean(cond4_all_elecs)  nanmean(cond2_all_elecs) nanmean(cond3_all_elecs)];
bar(bar_vectora1, 'b')
hold on
bar_vectora1sem =[nanstd(cond1_all_elecs)/sqrt(cond1_all_elecs_ntrls) nanstd(cond4_all_elecs)/sqrt(cond4_all_elecs_ntrls) nanstd(cond2_all_elecs)/sqrt(cond2_all_elecs_ntrls)...
    nanstd(cond3_all_elecs)/sqrt(cond3_all_elecs_ntrls) ];
errorbar(1:4,bar_vectora1,bar_vectora1sem, 'kx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'inseq +' 'outseq -' 'outseq +' 'inseq -' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])
title (band_value)
hgexport(gcf, [plot_type '_Barplots_SingleElec_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, [plot_type '_Barplots_SingleElec_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')

% reps  = 1000;
% adata = cond1_all_elecs;
% bdata = cond4_all_elecs;
% p1     = permutation_unpaired(adata, bdata, reps)
%
% adata = cond1_all_elecs;
% bdata = cond2_all_elecs;
% p2    = permutation_unpaired(adata, bdata, reps)
%
% adata = cond1_all_elecs;
% bdata = cond3_all_elecs;
% p3    = permutation_unpaired(adata, bdata, reps)
%
% reps  = 1000;
% adata = cond2_all_elecs;
% bdata = cond3_all_elecs;
% p4     = permutation_paired(adata, bdata, reps)
%
% adata = cond2_all_elecs;
% bdata = cond4_all_elecs;
% p5    = permutation_unpaired(adata, bdata, reps)
%
% adata = cond3_all_elecs;
% bdata = cond4_all_elecs;
% p6    = permutation_unpaired(adata, bdata, reps)
%
%
% alpha = .05
% pvals = [p1 p2 p3 p4 p5 p6];
% [p_fdr, p_masked] = fdr( pvals, alpha)
% disp(['inseq+ vs outseq- : ' num2str(pvals(3))])
% disp(['inseq+ vs outseq+ : ' num2str(pvals(1))])
% disp(['inseq+ vs inseq- : ' num2str(pvals(2))])
% disp(['outseq+ vs outseq- : ' num2str(pvals(5))])
% disp(['inseq- vs outseq- : ' num2str(pvals(6))])
% disp(['outseq+ vs inseq- : ' num2str(pvals(4))])
% p_fdr

mn = -.2
mx = .157
% bar plot for corr vs. incorr
corr_vect   = [cond1_all_elecs cond2_all_elecs];
incorr_vect = [cond3_all_elecs cond4_all_elecs];
corr_vect_ntrls   = length(corr_vect) - sum(isnan(corr_vect));
incorr_vect_ntrls = length(incorr_vect) - sum(isnan(incorr_vect));

figure
bar_vectora1 = [nanmean(corr_vect) nanmean(incorr_vect)];
bar(bar_vectora1, 'b')
hold on
bar_vectora1sem =[nanstd(cond1_all_elecs)/sqrt(cond1_all_elecs_ntrls) nanstd(cond2_all_elecs)/sqrt(cond2_all_elecs_ntrls)];
errorbar(1:2,bar_vectora1,bar_vectora1sem, 'kx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'Corr' 'Incorr' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])
title (band_value)
hgexport(gcf, [plot_type '_Barplots_CorrIncorr_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, [plot_type '_Barplots_CorrIncorr_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')

% bar plot for match vs. mismatch
match_vect    = [cond1_all_elecs cond4_all_elecs];
mismatch_vect = [cond2_all_elecs cond3_all_elecs];
match_vect_ntrls   = length(match_vect) - sum(isnan(match_vect));
mismatch_ntrls = length(mismatch_vect) - sum(isnan(mismatch_vect));

figure
bar_vectora1 = [nanmean(match_vect) nanmean(mismatch_vect)];
bar(bar_vectora1, 'b')
hold on
bar_vectora1sem =[nanstd(match_vect)/sqrt(match_vect_ntrls) nanstd(mismatch_vect)/sqrt(mismatch_ntrls)];
errorbar(1:2,bar_vectora1,bar_vectora1sem, 'kx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'match' 'mismatch' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])
title (band_value)
hgexport(gcf, [plot_type '_Barplots_MatchMismatch_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, [plot_type '_Barplots_MatchMismatch_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')


% bar plot for InSeq vs. OutSeq +
figure
bar_vectora1 = [nanmean(cond1_all_elecs) nanmean(cond2_all_elecs)];
bar(bar_vectora1, 'b')
hold on
bar_vectora1sem =[nanstd(corr_vect)/sqrt(corr_vect_ntrls) nanstd(incorr_vect)/sqrt(incorr_vect_ntrls)];
errorbar(1:2,bar_vectora1,bar_vectora1sem, 'kx')%ylim([mn mx])
set(gca, 'XTick', 1:4, 'XTickLabel', {'InSeq' 'OutSeq' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
clear bar_vectora1 bar_vectora1std
ylim([mn mx])
title (band_value)
hgexport(gcf, [plot_type '_Barplots_InSeqCorr_OutSeqCorr_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.eps'], hgexport('factorystyle'), 'Format', 'eps')
hgexport(gcf, [plot_type '_Barplots_InSeqCorr_OutSeqCorr_'  band_value '_' lock '_' num2str(start_time) '_to_' num2str(end_time) '.jpeg'], hgexport('factorystyle'), 'Format', 'jpeg')



