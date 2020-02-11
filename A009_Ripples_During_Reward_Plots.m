% plotting ripples within trials

for chan =1:chan_counter
    data = load([chan_name num2str(chan_length(chan)) '.mat']); % load exmp chan
    
    % extract time-locked LFP data
    LFP = ExtractTrialData_SM(pokeInAlignedBehavMatrix, filtfilt(flt,data.statMatrix(:,2)) );
    
    % extract trials' time-locked LFP data
    LFP_1 = cell2mat(LFP(trial_log_1))';
    
    figure
    for i = 1:5
        subplot(5,1,i)
        plot(LFP_1(i,400:900))
        xlim([0 length(LFP_1(i,400:900))])
    end
     title(['chan ' num2str(chan_length(chan))])
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
end