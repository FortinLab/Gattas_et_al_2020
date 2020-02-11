clear all; close all;clc
anim            = 6
task            = 'welltrained' %novel1, novel2, welltrained
fs              = 1000;
lock            = 'p_o'
%%
if anim == 1
    if strcmp('novel1',task)
        chan_name = 'SuperChris-Novel1-IntBadCutLFP_T';
    elseif strcmp('novel2',task)
        chan_name = 'SuperChris-Novel2-IntBadCut.plx_T';
    else
        chan_name = 'SuperChris-2-12-09_SG_final_T';
    end
    chan_length = [1:10 12:16 18:23];
elseif anim ==2
    if strcmp('novel1',task)
        chan_name = 'Stella-Novel1-IntBadCutLFP_T';
        chan_length = [2:10 12:23];
    elseif strcmp('novel2',task)
        chan_name = 'Stella-Novel2-IntBadCutLFP_T';
        chan_length = [2:10 12:23];
    else
        chan_name = 'Stella-2-12-2009_SG_final_T';
        chan_length = [2:10 12:23];
    end
elseif anim==3
    if strcmp('novel1',task)
        chan_name = 'Barat-Novel1-IntBadCutLFP_T';
    elseif strcmp('novel2',task)
        chan_name = 'Barat-Novel2-IntBadCutLFP_T';
    else
        chan_name = 'Barat-11-06-2008Skips_mrg_SG_final_T';
    end
    chan_length  = [1 3:10 12:21 23];
elseif anim==4
    if strcmp('novel1',task)
        chan_name = 'Buchanan-Novel1-IntBadCutLFP_T';
    elseif strcmp('novel2',task)
        chan_name = 'Buchanan-Novel2-IntBadCutLFP_T';
    else
        chan_name ='Buchanan4-20-withskips_mrg_SG_final_T';
    end
    chan_length = [1 2 4:10 12:13 15:23];
elseif anim ==5
    chan_name ='Mitt_July18_5odorswithSkips_SG_final_T';
    chan_length = [1:10 12:23];
elseif anim ==6
    chan_name = 'SAS01_SessiongGE44_mrg_GEcut1_T';
    chan_length = [1:2 4:10 12:18 21:24];
end
chan_counter = length(chan_length);
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])
load('baseline_wavelet_32num_20logdb_3hz_250hz_welltrained_avg_freq_ODOR_RUN_states')
load('BehaviorMatrix.mat')
load([chan_name '4.mat']) % load exmp chan

%% spatial map indications for each anim
[medial_elecs, lateral_elecs, i, i2,...
 subplot_row_lateral, subplot_row_medial , ...
 subplot_column_lateral, subplot_column_medial] = electrode_map( anim );

%%

% define eventWindow/lock
if strcmp('p_i', lock)
    eventWindow = [-0.5 1];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn');
    
elseif strcmp('p_o', lock)
    eventWindow = [-0.5 0.5];
    pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut');
end
inSeqLog     = [pokeInAlignedBehavMatrix.TranspositionDistance]==0;
corrTrlLog   = [pokeInAlignedBehavMatrix.Performance]==1;
odor_excld_A = [pokeInAlignedBehavMatrix.Odor]~=1;
inSeqCorrLog = inSeqLog&corrTrlLog&odor_excld_A==1;

%%
trial_log_1 = inSeqCorrLog;
LFP   = ExtractTrialData_SM(pokeInAlignedBehavMatrix, statMatrix(:,2));
LFP_1 = cell2mat(LFP(trial_log_1))';

%% create notch filter
fs=1000;
flt = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);
%% overlap elecs spectrum
% F           = 2.^(1.6:0.08:8);
% figure
% hold on
% for chan = 1:chan_counter
%      % LFP data
%     data = load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
%     
%     % extract time-locked LFP data
%     LFP = ExtractTrialData_SM(pokeInAlignedBehavMatrix, filtfilt(flt,data.statMatrix(:,2)) );
%     
%     % extract trials' time-locked LFP data
%     LFP_1 = cell2mat(LFP(trial_log_1));
%     % load data
%     [pxx,f] = pwelch(LFP_1, [],[],F, fs);
%     avg_pxx = mean(pxx,2);
%     plot(f, avg_pxx)
%     xlim([3 100])
%     clear LFP_1
% end
% % legend({'chans 1:5','chans 6:10, 18','chans 13:14 16 19:23'})
% xlabel('freq')
% ylabel('power')
% set(gca, 'FontSize',14,'FontWeight','bold')
% title('power distribution during odor presentation')
%%
i_counter = 0;
figure
F           = 2.^(1.6:0.08:8);

ylim_value = [4.2*10^-3  3*10^-3  2.6*10^-3 .8*10^-3  7*10^-4  5.5*10^-3]

for elec = lateral_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_lateral(anim), subplot_column_lateral(anim), i2(i_counter))
   % LFP data
    data = load([chan_name num2str(elec) '.mat']); % load exmp chan
    
    % extract time-locked LFP data
    LFP = ExtractTrialData_SM(pokeInAlignedBehavMatrix, filtfilt(flt,data.statMatrix(:,2)) );
    
    % extract trials' time-locked LFP data
    LFP_1 = cell2mat(LFP(trial_log_1));
    % load data
    [pxx,f] = pwelch(LFP_1, [],[],F, fs);
    avg_pxx = mean(pxx,2);
    plot(f, avg_pxx)
    ylim([0 ylim_value(anim)])
    xlim([4 57])
    xlabel('freq')
    ylabel('power')
    title(['chan ' num2str(elec)])
    set(gca, 'FontSize',10,'FontWeight','bold')
    clear pxx
end
% spatial plots 1/f: medial

i_counter = 0;
figure
for elec = medial_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_medial(anim), subplot_column_medial(anim), i(i_counter))
  % LFP data
    data = load([chan_name num2str(elec) '.mat']); % load exmp chan
    
    % extract time-locked LFP data
    LFP = ExtractTrialData_SM(pokeInAlignedBehavMatrix, filtfilt(flt,data.statMatrix(:,2)) );
    
    % extract trials' time-locked LFP data
    LFP_1 = cell2mat(LFP(trial_log_1));

    [pxx,f] = pwelch(LFP_1, [],[],F, fs);
    avg_pxx = mean(pxx,2);
    plot(f, avg_pxx)
    ylim([0 ylim_value(anim)])
    xlim([4 57])
    xlabel('freq')
    ylabel('power')
    title(['chan ' num2str(elec)])
    set(gca, 'FontSize',10,'FontWeight','bold')
    clear pxx
end
