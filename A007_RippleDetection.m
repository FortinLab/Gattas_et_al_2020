clear all; close all
anim = 1
lock = 'p_o' %'p_o' % p_i: poke-in % p_o: poke-out
match = 'yes';
plot_type = '';

task = 'welltrained' %novel1, novel2, welltrained
per = 0.5; % caxis variable for plotting
scale = 'db'
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

%% load behavioral and neural data
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load([chan_name '5.mat']) % load exmp chan
load('baseline_info_wavelet_32num_20logdb_3hz_250hz_notched_artifact_reject.mat')
load('chan_artifact_thresh.mat')
load('BehaviorMatrix.mat')
load([chan_name '4.mat']) % load exmp chan
%%
fs = 1000;
tvec =0:1/fs: 1/fs * (size(statMatrix(:,2),1)-1);

% artifact and detection params
DHPC_artifact_threshold_multiplier = 2.5;
DHPC_ripple_threshold_multiplier = 2;
win = round(0.008*fs); % window for calculating ripple rms
interripple_interval_min = 0.05*fs;
asd = 0.75; % fraction of thresh for spread of ripple in left direction
bsd = 0.75; % fraction of thresh for spread of ripple in right direction
minimum_ripple_duration_multiplier = 1;

before_event = 1; %1 secondgggg
after_event = 1;  %1 second

% Designing the filter %
bpFilt = designfilt('bandpassfir','FilterOrder',500, ...
    'CutoffFrequency1',150,'CutoffFrequency2',250, ...
    'SampleRate',1000);

%% init mtx
chan_ripple_indices = cell(length(chan_length),1);
chan_ripple_cntr    = cell(length(chan_length),1);
ripple_mtx          = cell(length(chan_length),1);
parfor chan = 1:chan_counter
    % LFP data
    data = load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
    LFP = data.statMatrix(:,2); % data
    
    % Check for LFP artifacts
    DHPC_artifact_threshold = [mean(LFP) + DHPC_artifact_threshold_multiplier*std(LFP),...
        mean(LFP) - DHPC_artifact_threshold_multiplier*std(LFP)];
    DHPC_artifacts_indx = find(LFP > DHPC_artifact_threshold(1) | LFP < DHPC_artifact_threshold(2));
    %figure; plot(0:1/fs:(length(LFP)-1)/fs,LFP,DHPC_artifacts_indx/fs,LFP(DHPC_artifacts_indx),'.r'); grid on; shg
    
    % use filter and calc rms
    DHPC_LFP_filtered_150_250 = filtfilt(bpFilt.Coefficients,1,LFP);
    DHPC_rms = conv(DHPC_LFP_filtered_150_250 .^2, ones(win,1)/win,'same'); % calculating the ripple rms
    DHPC_thresh = mean(DHPC_rms(~ismember(1:length(DHPC_rms),...
        [DHPC_artifacts_indx-fs/200; DHPC_artifacts_indx; DHPC_artifacts_indx+fs/200])))...
        + DHPC_ripple_threshold_multiplier*std(DHPC_rms(~ismember(1:length(DHPC_rms),...
        [DHPC_artifacts_indx-fs/200; DHPC_artifacts_indx; DHPC_artifacts_indx+fs/200])));
    
    lsd = DHPC_thresh*asd;
    rsd = DHPC_thresh*bsd;
    
    [p,ip] = findpeaks(DHPC_rms,'minpeakheight',DHPC_thresh);
    len = length(p);
    S = [];
    E = [];
    
    for ii = 1:len
        jj = ip(ii);
        while jj > 0 && DHPC_rms(jj) > lsd
            jj = jj - 1;
        end
        if jj
            S =union(S, jj);
        else
            S = union(S, ip(1));
        end
        jj = ip(ii);
        while jj < length(tvec) && DHPC_rms(jj) > rsd
            jj = jj + 1;
        end
        if jj < length(tvec)
            E = union(E, jj);
        else
            E = union(E, ip(end));
        end
    end
    
    DHPC_ripples_indx =[S', E'];
    
    minimum_ripple_duration = mean(DHPC_ripples_indx(:,2) - DHPC_ripples_indx(:,1)) + ...
        minimum_ripple_duration_multiplier*std(DHPC_ripples_indx(:,2) - DHPC_ripples_indx(:,1));
    DHPC_ripples_indx(DHPC_ripples_indx(:,2) - DHPC_ripples_indx(:,1) < minimum_ripple_duration,:) = [];
    ripples_midpoints = mean(DHPC_ripples_indx, 2);
    a = find((ripples_midpoints(2:end) - ripples_midpoints(1:end-1,:)) < interripple_interval_min);
    DHPC_ripples_indx(a,2) = DHPC_ripples_indx(a+1,2);
    DHPC_ripples_indx(a+1,:) = [];
    
    DHPC_ripples_centers_indx = zeros(size(DHPC_ripples_indx,1),1);
    for ii =1:size(DHPC_ripples_indx,1)
        [x, y] = findpeaks( -DHPC_LFP_filtered_150_250(DHPC_ripples_indx(ii,1):DHPC_ripples_indx(ii,2))); %detecting ripple troughs
        DHPC_ripples_centers_indx(ii) = y(x == max(x) )+ DHPC_ripples_indx(ii,1) - 1; %index of largest trough
    end
    
    non_empty = [];
    for ii= 1:length(DHPC_ripples_centers_indx)
        if DHPC_ripples_centers_indx(ii)-before_event*fs > 0 &&...
                DHPC_ripples_centers_indx(ii)+after_event*fs < length(LFP) && ...
                isempty(intersect(DHPC_ripples_centers_indx(ii)-before_event*fs:DHPC_ripples_centers_indx(ii)+after_event*fs, ...
                DHPC_artifacts_indx))
            non_empty = [non_empty ii];
        end
    end
    DHPC_ripples_centers_indx = DHPC_ripples_centers_indx(non_empty);
    DHPC_ripples_indx = DHPC_ripples_indx(non_empty,:);
    chan_ripple_indices{chan} = DHPC_ripples_indx;
    chan_ripple_cntr{chan} = DHPC_ripples_centers_indx;
    for a = 1:length(DHPC_ripples_centers_indx)
        ripple_mtx{chan}(a,:) = [LFP(abs(fs*1-DHPC_ripples_centers_indx(a)):abs(fs*1+DHPC_ripples_centers_indx(a)))];
    end
end

cd ./ripples
save('ripple_times', 'chan_ripple_indices', 'chan_ripple_cntr', 'ripple_mtx')
%% plot all detected ripples
for chan = 1:chan_counter
    % LFP data
    data = load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
    LFP = data.statMatrix(:,2); % data
    figure;
    ax1 = subplot(2,1,1);
    title(['thrshld_mul_' num2str(DHPC_ripple_threshold_multiplier)]);
    plot(0:1/fs:(length(LFP)-1)/fs, LFP); grid on; hold on;
    plot(DHPC_ripples_centers_indx/fs, LFP(DHPC_ripples_centers_indx),'ok');
    plot(DHPC_ripples_indx(:,1)/fs, LFP(DHPC_ripples_indx(:,1)),'*r',DHPC_ripples_indx(:,2)/fs, LFP(DHPC_ripples_indx(:,2)), '*c');
    ax2 = subplot(2,1,2);
    plot(0:1/fs:(length(LFP)-1)/fs, DHPC_rms); grid on; hold on;
    plot(DHPC_ripples_indx(:,1)/fs, DHPC_rms(DHPC_ripples_indx(:,1)),'*r',DHPC_ripples_indx(:,2)/fs, DHPC_rms(DHPC_ripples_indx(:,2)), '*c');
    plot(DHPC_ripples_centers_indx/fs, DHPC_rms(DHPC_ripples_centers_indx),'ok');
    line([0 (length(LFP)-1)/fs],[DHPC_thresh DHPC_thresh],'color', 'r');
    linkaxes([ax1 ax2], 'x');
end