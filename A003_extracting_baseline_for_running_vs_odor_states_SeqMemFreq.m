%% measure baseline
clear all; close all; clc
anim  = 6
task  = 'welltrained'
scale = 'db'
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load([chan_name '5.mat']) % load exmp chan
load('chan_artifact_thresh.mat')

% wavelet params
fs = 1000;
dt = 1/fs;
NumVoices = 32; %? %used 64 prev
a0 = 2^(1/NumVoices); %?
wavCenterFreq = 6/(2*pi); %?
minfreq = 3;
maxfreq = 250;
minscale = wavCenterFreq/(maxfreq*dt); %?
maxscale = wavCenterFreq/(minfreq*dt); %?
minscale = floor(NumVoices*log2(minscale));%?
maxscale = ceil(NumVoices*log2(maxscale)); %?
scales = a0.^(minscale:maxscale).*dt;
freq = wavCenterFreq./(fs*scales.*dt);

% filter
fs = 1000;
flt = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);

%
load([chan_name num2str(chan_length(2)) '.mat'])
baseline2 = cell(1,chan_counter);
data = single(zeros(chan_counter,length(statMatrix(:, 2))));
for chan = 1:chan_counter
    load([chan_name num2str(chan_length(chan)) '.mat']); % load chan
    data(chan,:) = single(filtfilt(flt,statMatrix(:, 2))); % notch filter
end
if size(statMatrix,1)>length(1:2000000);
    data = data(:, 1:2000000); %33 min
end

artifact_points = zeros(chan_counter, size(statMatrix,1));
chan_artifact_thresh = zeros(chan_counter,2);
for chan =1:chan_counter
    %load chan
    load([chan_name num2str(chan_length(chan))])
    
    % filter data
    x = filtfilt(flt,statMatrix(:, 2));
    
    %get mean and std
    mn = mean(x);
    sd = std(x,1);
    pos_thresh = mn+(5*sd);
    neg_thresh = mn-(5*sd);
    idx_above_thresh = (x>pos_thresh);
    idx_below_thresh = (x<neg_thresh);
 
    % get artifact indices/chan
    temp = idx_above_thresh+idx_below_thresh;
    artifact_idx = find(temp);
    for counter = 1:length(artifact_idx)
        if artifact_idx(counter)<1000
         temp(1:(artifact_idx(counter)+1000)) = 1; %set all points before and 1 sec after = artifact

        elseif artifact_idx(counter)>size(statMatrix,1)-1000
         temp((artifact_idx(counter)-1000):end) = 1; %set sec bef and sec after = artifact

        else
        temp((artifact_idx(counter)-1000):(artifact_idx(counter)+1000)) = 1; %set sec bef and sec after = artifact
        end
    end
    
    artifact_points(chan,:) = temp;  
    chan_artifact_thresh(chan,:) = [neg_thresh pos_thresh];
    clear temp
end

% extract freq info
parfor chan = 1:chan_counter
    cwt = cwtft({data(chan, :),dt},...
        'scales',scales,'wavelet','morl');
    clean_points = artifact_points(chan, 1:size(data,2)) == 0;
    baseline2 {chan}=  single(20*log10(abs(cwt.cfs(:,clean_points)).^2));
end


% 5 bands of interest
freq_idx (:,1) = freq>6 & freq<9;
freq_idx (:,2) = freq>9 & freq<12;
freq_idx (:,3) = freq>11 & freq<19;
freq_idx (:,4) = freq>19 & freq<36;
freq_idx (:,5) = freq>36 & freq<55;
freq_idx (:,6) = freq>60 & freq<100;
freq_idx (:,7) = freq>126 & freq<250;


% first average the freq comp of each band
freq_range_mean_chan = cell(1,chan_counter);
for elec = 1:chan_counter
    for f = 1:size(freq_idx,2) % loop thru the freq of that band
    freq_range_mean_chan{elec}(f,:) = mean(baseline2{elec}(freq_idx(:,f), :), 1);
    end
end

freq_std_chan = zeros(chan_counter, size(freq_idx,2));
freq_mean_chan = zeros(chan_counter, size(freq_idx,2));
for elec =1:chan_counter
    freq_std_chan(elec,:)  = std(freq_range_mean_chan{elec}, 1, 2);
    freq_mean_chan(elec,:) = mean(freq_range_mean_chan{elec},2);
end

save(['baseline_wavelet_32num_20logdb_3hz_250hz_' task '_avg_freq_ODOR_RUN_states_SeqMemFreq'], 'freq_std_chan', 'freq_mean_chan');