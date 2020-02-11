reps = 100000;
elec  = 1

adata = cond1(elec,:);
bdata = cond2(elec,:);
p1 = permutation_paired(adata, bdata, reps)

adata = cond1(elec,1:5);
bdata = cond3(elec,1:5);
p2 = permutation_paired(adata, bdata, reps)

adata = cond1(elec,:);
bdata = cond4(elec,:);
p3 = permutation_paired(adata, bdata, reps)

adata = cond2(elec,:);
bdata = cond3(elec,:);
p4 = permutation_paired(adata, bdata, reps)

adata = cond2(elec,:);
bdata = cond4(elec,:);
p5 = permutation_paired(adata, bdata, reps)

adata = cond3(elec,:);
bdata = cond4(elec,:);
p6 = permutation_paired(adata, bdata, reps)

%%
[p_fdr, p_masked] = fdr([p1 p2 p3 p4 p5 p6], 0.05)
                            % all poss p-values