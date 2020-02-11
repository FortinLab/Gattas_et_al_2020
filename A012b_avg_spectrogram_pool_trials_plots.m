% this scrip will 1_ plot running spectrogram across animals using pooled
% trias. No need to run other script to load variables
trial_selection = 'run'
pool_trials     = ''
% change dir to where indiv anim data is
if strcmp('run', trial_selection)
    cd(['D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim\fig2_odorVsrun\run_matrices_AL_PM'])
end
anim_list   = [1:5];
anim_length = length(anim_list);
% init
cond1_ch1_temp = cell(1,anim_length);
cond1_ch2_temp = cell(1,anim_length);
cond1_ch3_temp = cell(1,anim_length);
cond1_ch4_temp = cell(1,anim_length);

% get all trials per 4 chans
if strcmp('yes',pool_trials)
    for anim = 1:anim_length
        load(['anim' num2str(anim_list(anim)) '_AL_PM_run.mat'])
        cond1_ch1_temp{anim} = run_spectrogram_AL_PM(:,:,:,1);
        cond1_ch2_temp{anim} = run_spectrogram_AL_PM(:,:,:,2);
        cond1_ch3_temp{anim} = run_spectrogram_AL_PM(:,:,:,3);
        cond1_ch4_temp{anim} = run_spectrogram_AL_PM(:,:,:,4);
    end
    % conc all trials for each chan
    cond1_ch1_temp_temp = cat(3,cond1_ch1_temp{:});
    cond1_ch2_temp_temp = cat(3,cond1_ch2_temp{:});
    cond1_ch3_temp_temp = cat(3,cond1_ch3_temp{:});
    cond1_ch4_temp_temp = cat(3,cond1_ch4_temp{:});
    
else
    load(['anim' num2str(anim_list(1)) '_AL_PM_run.mat'])
    cond1_ch1_temp_temp = zeros(size(run_spectrogram_AL_PM,1), size(run_spectrogram_AL_PM,2), anim_length);
    cond1_ch2_temp_temp = zeros(size(run_spectrogram_AL_PM,1), size(run_spectrogram_AL_PM,2), anim_length);
    cond1_ch3_temp_temp = zeros(size(run_spectrogram_AL_PM,1), size(run_spectrogram_AL_PM,2), anim_length);
    cond1_ch4_temp_temp = zeros(size(run_spectrogram_AL_PM,1), size(run_spectrogram_AL_PM,2), anim_length);
    for anim = 1:anim_length
        load(['anim' num2str(anim_list(anim)) '_AL_PM_run.mat'])
        cond1_ch1_temp_temp(:,:,anim) = nanmean(run_spectrogram_AL_PM(:,:,:,1),3);
        cond1_ch2_temp_temp(:,:,anim) = nanmean(run_spectrogram_AL_PM(:,:,:,2),3);
        cond1_ch3_temp_temp(:,:,anim) = nanmean(run_spectrogram_AL_PM(:,:,:,3),3);
        cond1_ch4_temp_temp(:,:,anim) = nanmean(run_spectrogram_AL_PM(:,:,:,4),3);
    end
    
end
% wavelet params
fs = 1000;
NumVoices = 32;
dt = 1/fs;
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

%% plot vars
edge = 100;
tickmarks = 1:20:length(freq);
mx = .8
mn = -.8

figure;
for i = 1:4
    if i ==1
        mtx = nanmean(cond1_ch1_temp_temp(:,edge+1:end-edge,:),3);
    elseif i ==2
        mtx = nanmean(cond1_ch2_temp_temp(:,edge+1:end-edge,:),3);
    elseif i ==3
        mtx = nanmean(cond1_ch3_temp_temp(:,edge+1:end-edge,:),3);
    elseif i ==4
        mtx = nanmean(cond1_ch4_temp_temp(:,edge+1:end-edge,:),3);
    end
    subplot(1,4,i)
    imagesc(linspace(-.5+edge/fs,.5-edge/fs,size(cond1_ch1_temp_temp,2)), 1:length(freq),  mtx)
    colorbar
    colormap jet
    caxis([-.75 1])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),  'FontSize', 16, 'FontWeight', 'bold')
    clear mtx
end

%% save figures all 4 elecs in 1 plot
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim\fig2_runVsOdor')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim\fig2_odorVsrun')
end
hgexport(gcf, ['group_AL_PM_' trial_selection '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
hgexport(gcf, ['group_AL_PM_' trial_selection '.eps'], hgexport('factorystyle'), 'Format', 'eps')
%% save figure - 1 elec at a time
for i =4
    figure;

    if i ==1
        mtx = nanmean(cond1_ch1_temp_temp(:,edge+1:end-edge,:),3);
    elseif i ==2
        mtx = nanmean(cond1_ch2_temp_temp(:,edge+1:end-edge,:),3);
    elseif i ==3
        mtx = nanmean(cond1_ch3_temp_temp(:,edge+1:end-edge,:),3);
    elseif i ==4
        mtx = nanmean(cond1_ch4_temp_temp(:,edge+1:end-edge,:),3);
    end
    imagesc(linspace(-.5+edge/fs,.5-edge/fs,size(cond1_ch1_temp_temp,2)), 1:length(freq),  mtx)
    colorbar
    colormap jet
    caxis([-.75 .75])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),  'FontSize', 16, 'FontWeight', 'bold')
    hgexport(gcf, ['group_AL_PM_' trial_selection '_elec_' num2str(i)  '.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
    hgexport(gcf, ['group_AL_PM_' trial_selection '_elec_' num2str(i)  '.eps'], hgexport('factorystyle'), 'Format', 'eps')

    clear mtx
end


%% make bar plots of gamma power


cond_1_mean = squeeze(nanmean(squeeze(nanmean(cond1_ch1_temp_temp(freq_range, :,:),1)),1));
cond_1_ntrls = length(cond_1_mean) - sum(isnan(cond_1_mean));

cond_2_mean = squeeze(nanmean(squeeze(nanmean(cond1_ch2_temp_temp(freq_range, :,:),1)),1));
cond_2_ntrls = length(cond_2_mean) - sum(isnan(cond_2_mean));

cond_3_mean = squeeze(nanmean(squeeze(nanmean(cond1_ch3_temp_temp(freq_range, :,:),1)),1));
cond_3_ntrls = length(cond_3_mean) - sum(isnan(cond_3_mean));

cond_4_mean = squeeze(nanmean(squeeze(nanmean(cond1_ch4_temp_temp(freq_range, :,:),1)),1));
cond_4_ntrls = length(cond_4_mean) - sum(isnan(cond_4_mean));

bar_vector1 = [nanmean(cond_1_mean) nanmean(cond_2_mean) nanmean(cond_3_mean) nanmean(cond_4_mean)];
figure;bar(bar_vector1, 'k')
hold on
bar_vectora1sem = [nanstd(cond_1_mean)/sqrt(cond_1_ntrls) nanstd(cond_2_mean)/sqrt(cond_2_ntrls) nanstd(cond_3_mean)/sqrt(cond_3_ntrls) nanstd(cond_4_mean)/sqrt(cond_4_ntrls)  ];
errorbar(1:4,bar_vector1,bar_vectora1sem, 'kx')%ylim([mn mx])

hgexport(gcf, ['barplots_gamma.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
hgexport(gcf, ['barplots_gamma.eps'], hgexport('factorystyle'), 'Format', 'eps')
%%
close all;clc;
range = 'hifreq' % theta hifreq
cd(['D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim\fig2_odorVsrun\run_matrices_AL_PM'])
if strcmp('theta', range)   
colorList    = {'r-o', 'c--o', 'g-o', 'y-o', 'm--o',};
colorListSEM = {'r', 'c--', 'g', 'y', 'm--'};
freq_range = freq>4 & freq<12;
elseif strcmp('hifreq', range)
colorList    = {'r-o', 'c--o', 'g-o', 'y-o', 'm-o',};
colorListSEM = {'r', 'c', 'g', 'y', 'm'};
freq_range = freq>24
end
% add average trace
chan1 = squeeze(nanmean(squeeze(nanmean(cond1_ch1_temp_temp(freq_range, :,:),1)),1));
chan2 = squeeze(nanmean(squeeze(nanmean(cond1_ch2_temp_temp(freq_range, :,:),1)),1));
chan3 = squeeze(nanmean(squeeze(nanmean(cond1_ch3_temp_temp(freq_range, :,:),1)),1));
chan4 = squeeze(nanmean(squeeze(nanmean(cond1_ch4_temp_temp(freq_range, :,:),1)),1));

figure;hold on
bar_vectora1 = [nanmean(chan1) nanmean(chan2) nanmean(chan3) nanmean(chan4)];
bar_vectora1sem = [nanstd(chan1)/sqrt(length(chan1)) nanstd(chan2)/sqrt(length(chan2))...
    nanstd(chan3)/sqrt(length(chan3)) nanstd(chan4)/sqrt(length(chan4))];
b1= bar(bar_vectora1,'k')
errorbar(1:4,bar_vectora1,bar_vectora1sem,'k*')


colors       = {'r', 'c', 'g', 'y', 'm'};
band_value = '>24 hz';
mn = -.8;
mx = .5;
for anim = 1:5
    load(['anim' num2str(anim_list(anim)) '_AL_PM_run.mat'])

    cond_1_mean = [];
   
    for elec = 1:4
        cond_1_mean(elec,:) = squeeze(nanmean(squeeze(nanmean(run_spectrogram_AL_PM(freq_range, :,:,elec),1)),1));
    end
    
     % traces for indiv animals
    plot(nanmean(cond_1_mean,2) , colorList{anim},'MarkerFaceColor', colors{anim})
    bar_vectora1 = nanmean(cond_1_mean,2);
    bar_vectora1sem = [nanstd(cond_1_mean,0,2)/sqrt(size(cond_1_mean,2)) ];
    errorbar(1:4,bar_vectora1,bar_vectora1sem,colorListSEM{anim})  
    
%   p = anova1(cond_1_mean')
%   reps  = 1000;
%   adata = cond_1_mean(1,:);
%   bdata = cond_1_mean(2,:);
%   p1     = permutation_unpaired(adata, bdata, reps)
%   
%   adata = cond_1_mean(1,:);
%   bdata = cond_1_mean(3,:);
%   p2    = permutation_unpaired(adata, bdata, reps)
%   
%   adata = cond_1_mean(1,:);
%   bdata = cond_1_mean(4,:);
%   p3    = permutation_unpaired(adata, bdata, reps)
%   
%   adata = cond_1_mean(2,:);
%   bdata = cond_1_mean(3,:);
%   p4    = permutation_unpaired(adata, bdata, reps)
%   
%   adata = cond_1_mean(2,:);
%   bdata = cond_1_mean(4,:);
%   p5    = permutation_unpaired(adata, bdata, reps)
%   
%   adata = cond_1_mean(3,:);
%   bdata = cond_1_mean(4,:);
%   p6    = permutation_unpaired(adata, bdata, reps)
%     
%   alpha = .05
%   pvals = [p1 p2 p3 p4 p5 p6];
%   [p_fdr, p_masked] = fdr( pvals, alpha)
end
ylim([mn mx])
if strcmp('yes',pool_trials)
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_material_pool_trials_and_indiv_anim\fig2_runVsOdor')
else
    cd('D:\Gattas\ephys_data_final\group_plots\paper_figure_single_value_per_anim\fig2_odorVsrun')
end
hgexport(gcf, ['gamma_indiv_anim_traces.jpg'], hgexport('factorystyle'), 'Format', 'jpeg')
hgexport(gcf, ['gamma_indiv_anim_traces.eps'], hgexport('factorystyle'), 'Format', 'eps')