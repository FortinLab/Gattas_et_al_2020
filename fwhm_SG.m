function [fullwidth_at_halfmax,I_peak, half_max_idx_a, half_max_idx_b] = fwhm_SG(sig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


[p,ix] = findpeaks(sig);
I_peak = ix(p == max(p));
pk = sig(I_peak);
half_max_idx_a = [];
half_max_idx_b = [];
fullwidth_at_halfmax = [];
% if peak is positive
if ~isempty(I_peak) && pk>0
    jj = I_peak;
    while jj > 1 && (sig(jj) - 0.5*max(p)) > 0
        jj = jj-1;
    end
    
    kk = I_peak;
    while kk < length(sig) && (sig(kk) -0.5*max(p)) > 0
        kk = kk+1;
    end
    
    if jj == 0
        jj = 1;
    elseif kk == length(sig)
        kk = length(sig);
    end
    
    I_half = [jj, kk];
    
    if jj == 1 || kk == length(sig) || max(p) <= 0
        fullwidth_at_halfmax = nan;
        
    else
        fullwidth_at_halfmax = mean(sig(jj:kk));
        half_max_idx_a = jj;
        half_max_idx_b = kk;
    end
    
    
    % if peak is negative
elseif ~isempty(I_peak) && pk<0
    % find the end point
    diff_val = diff(sig(I_peak:end));
    % find first 4 consec positive points
    pos_val = find(diff_val>0);
    dif_pos_val = diff(pos_val);
    for a = 1:length(dif_pos_val)
        if sum(dif_pos_val(a:a+3))==4 % end point is before 4 poitns of positive deflection
            baseline = pos_val(a);
            return
        end
    end
    end_point = I_peak+baseline;
    
    % find the begining point
    diff_val = diff(sig(1:I_peak));
    % find first 4 consec neg points
    neg_val = find(diff_val<0);
    dif_neg_val = fliplr(diff(neg_val));
    for a = 1:length(dif_neg_val)
        if sum(dif_neg_val(a:a+3))==4 % end point is before 4 poitns of positive deflection
            baseline = a;
            return
        end
    end
    begin_point =neg_val(end-baseline);
    
    % find half the peak value
    half_peak = 0.5*(abs((sig(end_point)+sig(begin_point))/2) - abs(sig(I_peak))) ;
    half_peak_value =sig(I_peak)- half_peak;
    
    % find half max to the left
    decay = sig(I_peak:end_point);
    idx = find(decay-half_peak_value>0);
    min_val =  min(decay(idx)-half_peak_value);
    
    for a = I_peak:end_point
        if sig(a)-half_peak_value == min_val
            half_max_idx_b = a;
            return
        end
    end
    
    % find half max to the right
    rise = sig(begin_point:I_peak);
    idx = find(rise - sig(half_max_idx_b)>0)
    min_val = min(rise(idx)-sig(half_max_idx_b));
    
    for a = begin_point:I_peak
        if sig(a)-sig(half_max_idx_b) == min_val
            half_max_idx_a = a;
            return
        end
    end
    fullwidth_at_halfmax = mean(sig(half_max_idx_a:half_max_idx_b));
    
else
    fullwidth_at_halfmax = nan;
end

end

