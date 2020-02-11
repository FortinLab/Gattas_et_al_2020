function [ pvals, p_fdr, p_masked] = TrialTypeStats( cond_1_mean, cond_2_mean, cond_3_mean, cond_4_mean )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

reps  = 1000;
alpha = .05;

adata = cond_1_mean;
bdata = cond_4_mean;
p1     = permutation_unpaired(adata, bdata, reps);

adata = cond_1_mean;
bdata = cond_2_mean;
p2    = permutation_unpaired(adata, bdata, reps);

adata = cond_1_mean;
bdata = cond_3_mean;
p3    = permutation_unpaired(adata, bdata, reps);

adata = cond_2_mean;
bdata = cond_3_mean;
p4     = permutation_paired(adata, bdata, reps);

adata = cond_2_mean;
bdata = cond_4_mean;
p5    = permutation_unpaired(adata, bdata, reps);

adata = cond_3_mean;
bdata = cond_4_mean;
p6    = permutation_unpaired(adata, bdata, reps);


pvals = [p1 p2 p3 p4 p5 p6];
[p_fdr, p_masked] = fdr( pvals, alpha)
disp(['inseq+ vs outseq- : ' num2str(pvals(1))])
disp(['inseq+ vs outseq+ : ' num2str(pvals(2))])
disp(['inseq+ vs inseq- : ' num2str(pvals(3))])
disp(['outseq+ vs inseq- : ' num2str(pvals(4))])
disp(['outseq+ vs outseq- : ' num2str(pvals(5))])
disp(['inseq- vs outseq- : ' num2str(pvals(6))])



end

