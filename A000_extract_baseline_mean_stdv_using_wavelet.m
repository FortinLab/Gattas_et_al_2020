%% measure baseline
clear all; close all; clc
for anim  =1:5
clearvars -except anim
run   = ''
task  = 'welltrained'
scale = 'db'
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])

%%
[chan_name, chan_length ]= get_anim_info(anim, task);
chan_counter = length(chan_length);
load([chan_name '5.mat']) % load exmp chan

%% wavelet params
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

%% artifact rejection
fs = 1000;
flt1 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);
flt2 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
    'DesignMethod','butter','SampleRate',fs);

load([chan_name num2str(chan_length(4)) '.mat'])
baseline2 = cell(1,chan_counter);
data = single(zeros(chan_counter,length(statMatrix(:, 2))));
for chan = 1:chan_counter
    load([chan_name num2str(chan_length(chan)) '.mat']); % load chan
    data(chan,:) = single(filtfilt(flt1,statMatrix(:, 2)));
    clear  temp
end

if strcmp('yes',run)
    data = data(:, 1:1200000);
elseif size(statMatrix,1)>length(1:2000000);
    data = data(:, 1:2000000); %33 min
end
%% get artifact logical vector and save channel artifact threshold
load([chan_name num2str(chan_length(chan)) '.mat']); % load chan
artifact_points = zeros(chan_counter, size(statMatrix,1));
chan_artifact_thresh = zeros(chan_counter,2);
for chan =1:chan_counter
    %load chan
    load([chan_name num2str(chan_length(chan))])
    
    % filter data
    x= filtfilt(flt1,statMatrix(:, 2));    

    %get mean and std
    mn = mean(x);
    sd = std(x,1);
    pos_thresh = mn+(5*sd);
    neg_thresh = mn-(5*sd);
    idx_above_thresh = (x>pos_thresh);
    idx_below_thresh = (x<neg_thresh);
    figure; plot(x)
    line([0 length(x)] ,[pos_thresh pos_thresh], 'Color','k','LineWidth',4)
    line([0 length(x)] ,[neg_thresh neg_thresh], 'Color','k','LineWidth',4)
    
    % get artifact indices/chan
    temp = idx_above_thresh+idx_below_thresh;
    artifact_idx = find(temp);
    for counter = 1:length(artifact_idx)
        if artifact_idx(counter)<1000
         temp(1:(artifact_idx(counter)+1000)) = 1; %set sec bef and sec after = artifact
        elseif artifact_idx(counter)+1000>length(temp)
            temp(artifact_idx(counter):end) = 1;
        else
        temp((artifact_idx(counter)-1000):(artifact_idx(counter)+1000)) = 1; %set sec bef and sec after = artifact
        end
    end
    
    artifact_points(chan,:) = temp';  
    chan_artifact_thresh(chan,:) = [neg_thresh pos_thresh];
    clear temp
end

save('chan_artifact_thresh', 'chan_artifact_thresh')
%% extract freq info
parfor chan = 1:chan_counter
    cwt = cwtft({data(chan, :),dt},...
        'scales',scales,'wavelet','morl');
    clean_points = artifact_points(chan, 1:size(data,2)) == 0;
    baseline2 {chan}=  single(10*log10(abs(cwt.cfs(:,clean_points)).^2));

end

%% calc mean and stdv
chan_powr_mn  = zeros(length(freq),chan_counter);
chan_powr_std = zeros(length(freq),chan_counter);
for chan =1:chan_counter
    chan_powr_mn(:, chan)  = nanmean(baseline2{chan}, 2);
    chan_powr_std(:, chan) = nanstd(baseline2{chan}, 1, 2);
end
save(['baseline_info_wavelet_' num2str(NumVoices) 'num_10logdb_3hz_250hz_notched_artifact_reject'], 'chan_powr_mn', 'chan_powr_std');

end

