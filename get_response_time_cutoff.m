function [ crossing ] = get_response_time_cutoff(anim, task, times,inSeqLog ,otSeqLog, trial_num )

% provide session name and reponse time vector, and it will provide a time cut
% off based on the intersection of the two distributions (welltrained), or
% at the 2.5thpercentile for novel1/2 respectively of all response
% times. 

if strcmp('novel1',task) || strcmp('novel2',task)

    crossing = prctile (times(inSeqLog),5);   
    
    figure;subplot(211);hold on;  h = histfit (times); h(1).FaceColor = 'm'; h(2).Color = 'k';
    line([crossing crossing], [0 length(times)], 'color', 'k', 'LineWidth', 2)
    subplot(212); plot(trial_num(inSeqLog), times(inSeqLog), 'r*') % look at response distribution
    hold on
    plot(trial_num(otSeqLog), times(otSeqLog), 'b*')
    line([0 length(trial_num)], [crossing crossing], 'color', 'k', 'LineWidth', 2)
    suptitle(['anim' num2str(anim) ', time thresh ' num2str(crossing)])
    

elseif strcmp('welltrained',task)
%figure
h = histfit (times(inSeqLog))
h(1).FaceColor = 'm';
h(2).Color = 'k';

hold on
I = histfit (times(otSeqLog))
I(1).FaceColor = 'b';
I(2).Color = 'k';

pd  = fitdist(times(inSeqLog)','Normal')
mu1 = pd.mu;
sigma1 = pd.sigma;

pd  = fitdist(times(otSeqLog)','Normal')
mu2 = pd.mu;
sigma2 = pd.sigma;

% save these from histfit function
load('xy_responsetimeInseqPDF.mat')
x1 = x;
load('xy_responsetimeOutseqPDF.mat')
x2 = x;

x_tot = sort([x1 x2])
dist1 =1/sqrt(2*pi*sigma1^2) * exp(-(x_tot-mu1).^2/2/sigma1^2);
dist2 =1/sqrt(2*pi*sigma2^2) * exp(-(x_tot-mu2).^2/2/sigma2^2);
figure;hold on; plot(x_tot,dist1); plot(x_tot,dist2)

% determine indices
str_idx = find(x_tot>x1(1)); str_idx = str_idx(1); % begining of inseq dist
end_idx = find(x_tot<mu1); end_idx = max(end_idx);

counter_vector = str_idx:end_idx;
figure;hold on; plot(x_tot,dist1); plot(x_tot,dist2);
plot(x_tot(counter_vector),dist1(counter_vector), 'k'); plot(x_tot(counter_vector),dist2(counter_vector), 'k')

diff = zeros(length(counter_vector),1);
for i = 1:length(counter_vector)
    i
    diff(i) = abs(dist1(counter_vector(i)) - dist2(counter_vector(i)));
end
crossing = x_tot(counter_vector(find(diff==min(diff))));


close all;
figure;subplot(211);hold on; plot(x_tot,dist1); plot(x_tot,dist2);
plot(x_tot(counter_vector),dist1(counter_vector), 'k'); plot(x_tot(counter_vector),dist2(counter_vector), 'k')
line([crossing crossing], [0 max([dist1 dist2])], 'color', 'k', 'LineWidth', 2)

subplot(212); plot(trial_num(inSeqLog), times(inSeqLog), 'r*') % look at response distribution
hold on
plot(trial_num(otSeqLog), times(otSeqLog), 'b*')
line([0 length(trial_num)], [crossing crossing], 'color', 'k', 'LineWidth', 2)
suptitle(['anim' num2str(anim) ', time thresh ' num2str(crossing)])


end



end

