
clear elec_4 elec_3 elec_2 elec_1 elec_power
elec_1(:,2) = [cond4_trials{1,1:6}];
elec_1(:,1) = ones;
elec_2(:,2) = [cond4_trials{2,1:6}];
elec_2(:,1) = 2;
elec_3(:,2) = [cond4_trials{3,1:6}];
elec_3(:,1) = 3;
elec_4(:,2) = [cond4_trials{4,1:6}];
elec_4(:,1) = 4;
elec_power = [elec_1;elec_2;elec_3;elec_4]
[r,p] = corrcoef(elec_power(:,1), elec_power(:,2))
corr_val = r(2);
sig = p(2)
[corr_val sig]
figure; plot(elec_power(:,1), elec_power(:,2), '*r')