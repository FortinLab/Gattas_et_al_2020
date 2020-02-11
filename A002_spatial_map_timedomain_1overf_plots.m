% spatial maps: time domain and 1/f plots
clear all; close all
anim = 4
task = 'welltrained' %novel1, novel2, welltrained
cd(['D:\Gattas\ephys_data_final\' task '\anim' num2str(anim)])
%% load behavioral and neural data
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
    elseif strcmp('novel2',task)
        chan_name = 'Stella-Novel2-IntBadCutLFP_T';
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
load('chan_artifact_thresh.mat')

%% spatial map indications for each anim
[medial_elecs, lateral_elecs, i, i2,...
 subplot_row_lateral, subplot_row_medial , ...
 subplot_column_lateral, subplot_column_medial] = electrode_map( anim );

%% filter, address artifacts, and save data in mtx

% create filter
fs = 1000;
flt = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);

% filter data and save
load([chan_name '4.mat']) % load exmp chan
data = (zeros(chan_counter,length(statMatrix(:, 2))));
for chan = 1:chan_counter
    % get data
    load([chan_name num2str(chan_length(chan)) '.mat']);
    
    % notch filter
    x = (filtfilt(flt,statMatrix(:, 2)));
    
    data(chan,:) = x;
end
data = single(data);
%% replace artifacts with interpolation

data_artifact_free = zeros(chan_counter,length(statMatrix(2:end-1, 2)));
for chan = 1:chan_counter
    x =  data(chan,:);
    % find idx above and bel threshold
    idx_above_thresh = x>chan_artifact_thresh(chan,2);
    idx_below_thresh = x<chan_artifact_thresh(chan,1);
    
    % total idx
    temp = idx_above_thresh+idx_below_thresh;
    artifact_idx = find(temp);
    end_pnt = size(data,2); %size of time series - alst point
    % loop thru artifacts and replace the w/ an avg of point before and
    % point after
    artifact_idx_value = zeros(length(artifact_idx),1);
    parfor counter = 1:length(artifact_idx)
        if artifact_idx(counter) == 1 || artifact_idx(counter) == end_pnt
            artifact_idx_value(counter) = nan;
            
        else
            % interporalte by averaging nearby points
            artifact_idx_value(counter) = mean([x((artifact_idx(counter)-1))  x((artifact_idx(counter)-1))]);
        end
        
    end
    data(chan,artifact_idx) = artifact_idx_value';
    
    clear artifact_idx artifact_idx_value
end
data = data(:,2:end-1);
%% plot all elecs overlaped
% F           = 2.^(1.6:0.08:8);
% figure
% hold on
% for elec = 1:length(chan_length)
%     % load data
%     [pxx,f] = pwelch(data(elec,:), [],[],F, fs);
%     plot(f, pxx, 'LineWidth', 1)
%     xlim([3 57])
% end
% xlabel('freq')
% ylabel('power')
% set(gca, 'FontSize',14,'FontWeight','bold')

%% group 1 elecs: superchris
F           = 2.^(1.6:0.08:8);

if anim==1
    figure
    hold on
    anim1_1 = [1:5];
    anim1_2 = [6:10 18];
    anim1_3 = [12 13:14 16 19:23];
    
    for elec = 1:length(anim1_1)
        % load data
        [pxx,f] = pwelch(data(find(anim1_1(elec)==chan_length),:), [],[],F, fs);
        plot(f, pxx, 'm', 'LineWidth',1)
        xlim([3 57])
    end
    
    % group 2 elecs
    for elec = 1:length(anim1_2)
        % load data
        [pxx,f] = pwelch(data(find(anim1_2(elec)==chan_length),:), [],[],F, fs);
        plot(f, pxx, 'k', 'LineWidth', 1)
        xlim([3 57])
    end
    
    % group 3 elecs
    for elec = 1:length(anim1_3)
        % load data
        [pxx,f] = pwelch(data(find(anim1_3(elec)==chan_length),:), [],[],F, fs);
        plot(f, pxx, 'b', 'LineWidth', 1)
        xlim([3 57])
    end
    % legend({'chans 1:5','chans 6:10, 18','chans 13:14 16 19:23'})
    xlabel('freq')
    ylabel('power')
    set(gca, 'FontSize',14,'FontWeight','bold')
    
elseif anim ==6
    anim6_1 = [1 2 5 7 8 10:14 16:19];
    anim6_2 = [3 4 6 9 15  20];
    
    figure
    hold on
    for elec = anim6_1
        % load data
        [pxx,f] = pwelch(data(elec,:), [],[],F, fs);
        plot(f, pxx, 'k', 'LineWidth',1)
        xlim([3 57])
    end
    
    % group 2 elecs
    for elec = anim6_2
        % load data
        [pxx,f] = pwelch(data(elec,:), [],[],F, fs);
        plot(f, pxx, 'r', 'LineWidth', 1)
        xlim([3 57])
    end
    xlabel('freq')
    ylabel('power')
    set(gca, 'FontSize',14,'FontWeight','bold')
    
end
%% spatial plots 1/f: lateral
F           = 2.^(1.6:0.08:8);
ylim_value  = [4.5*10^-3 3.5*10^-3 2.5*10^-3 3.4*10^-3 6.8*10^-3 .03] 
i_counter = 0;
figure
for elec = lateral_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_lateral(anim), subplot_column_lateral(anim), i2(i_counter))
    [pxx,f] = pwelch(data(find(elec==chan_length),:), [],[],F, fs);
    plot(f, pxx, 'm', 'LineWidth',4)
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
    [pxx,f] = pwelch(data(find(elec==chan_length),:), [],[],F, fs);
    plot(f, pxx, 'm', 'LineWidth',4)
    ylim([0 ylim_value(anim)])
    xlim([4 57])
    xlabel('freq')
    ylabel('power')
    title(['chan ' num2str(elec)])
    set(gca, 'FontSize',10,'FontWeight','bold')
    clear pxx
end

%% spatial plots time domain data: lateral and medial
i_counter = 0;
figure
for elec = lateral_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_lateral(anim), subplot_column_lateral(anim), i2(i_counter))
    plot(data(find(elec==chan_length),9000:11000), 'k')
    ylim([-1 1])
    xlim([0 2000])
    xlabel('time')
    ylabel('power')
    title(['chan ' num2str(elec)])
    set(gca, 'FontSize',10,'FontWeight','bold')
end

i_counter = 0;
figure
for elec = medial_elecs
    i_counter = i_counter+1;
    subplot(subplot_row_medial(anim), subplot_column_medial(anim), i(i_counter))
    plot(data(find(elec==chan_length),9000:11000), 'k')
    ylim([-1 1])
    xlim([0 2000])
    xlabel('time')
    ylabel('power')
    title(['chan ' num2str(elec)])
    set(gca, 'FontSize',10,'FontWeight','bold')
end