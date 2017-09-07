function submit_parafac_sub10(dd, nclass, nsub, nrun, nburn, nthin)

cvidx = repmat(1:10, nsub, 1);
cvidx = cvidx(:);

subidx = repmat(1:nsub, 1, 10);
subidx = subidx(:);

nrepf = cvidx(dd);
nsubf = subidx(dd);

disp(['nrep: ', num2str(nrepf) ' nsub: ', num2str(nsubf)]);

load(strcat('/Shared/ssrivastva/wasp/parafac/data/sub10/parafac_cv_', num2str(nrepf), '_sub_', num2str(nsubf), '_k10.mat'));

train = partMat;
[nsample ndim] = size(train);
cats = repmat(2, 1, ndim);

[history, tend] = parafac_dx_sub(train, cats, nclass, 100000, nrun, nburn, nthin);

save(strcat('/Shared/ssrivastva/wasp/parafac/result/sub10/samp/res_cv_', num2str(nrepf), '_sub_', num2str(nsubf), '_k10.mat'), 'history', 'tend');    
csvwrite(strcat('/Shared/ssrivastva/wasp/parafac/result/sub10/samp/time_cv_', num2str(nrepf), '_sub_', num2str(nsubf), '_k10.csv'), tend);

disp(['done with cv ' num2str(nrepf) ' ...' ' subset ... ' num2str(nsubf) ' ... ']);    
        
quit


