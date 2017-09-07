function submit_parafac_full(dd, nclass, nrun, nburn, nthin)

load('/Shared/ssrivastva/wasp/parafac/data/parafac_full_data.mat');

train = Yn_big(:, :, dd); 
[nsample ndim] = size(train);
cats = repmat(2, 1, ndim);
[history, tend] = parafac_dx_com(train, cats, nclass, nrun, nburn, nthin);

save(strcat('/Shared/ssrivastva/wasp/parafac/result/full/res_', num2str(dd), '.mat'), 'history', 'tend');    
csvwrite(strcat('/Shared/ssrivastva/wasp/parafac/result/full/time_', num2str(dd), '.csv'), tend);

disp(['done with rep ' num2str(dd) ' ...' ]);    
        
quit
