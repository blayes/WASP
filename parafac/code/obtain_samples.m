% mcmc samples
clear; 
for cc=1:10
    load(strcat('/Shared/ssrivastva/wasp/parafac/result/full/res_', num2str(cc), ...
                '.mat'));    
    for dd = 1:20
        margMat = zeros(1000, 2);
        for ss = 1:1000
            margMat(ss, :) = history{1, ss}(dd, :);
        end
        csvwrite(strcat('/Shared/ssrivastva/wasp/parafac/result/full/res_cv_', num2str(cc), ...
            '_dim_', num2str(dd), '.csv'), margMat);    
    end   
end

% wasp samples 
clear;
for cc = 1:10
    for kk=1:5
        load(strcat('/Shared/ssrivastva/wasp/parafac/result/sub5/samp/res_cv_', ...
                    num2str(cc), '_sub_', num2str(kk), '_k5.mat'));    
        for dd = 1:20
            margMat = zeros(1000, 2);
            for ss = 1:1000
                margMat(ss, :) = history{1, ss}(dd, :);
            end
            csvwrite(strcat('/Shared/ssrivastva/wasp/parafac/result/sub5/samp/csv/res_cv_', num2str(cc), ...
                            '_sub_', num2str(kk), '_dim_', num2str(dd), '_k5.csv'), margMat);    
        end   
    end    
end

clear;
for cc = 1:10
    for kk=1:10
        load(strcat('/Shared/ssrivastva/wasp/parafac/result/sub10/samp/res_cv_', ...
                    num2str(cc), '_sub_', num2str(kk), '_k10.mat'));    
        for dd = 1:20
            margMat = zeros(1000, 2);
            for ss = 1:1000
                margMat(ss, :) = history{1, ss}(dd, :);
            end
            csvwrite(strcat('/Shared/ssrivastva/wasp/parafac/result/sub10/samp/csv/res_cv_', num2str(cc), ...
                            '_sub_', num2str(kk), '_dim_', num2str(dd), '_k10.csv'), margMat);    
        end   
    end    
end

% cmc; sdp samples 

clear; 
for cc = 1:10
    for kk=1:5
        load(strcat('/Shared/ssrivastva/wasp/parafac/result/comp/sub5/samp/res_cv_', ...
                    num2str(cc), '_sub_', num2str(kk), '_k5.mat'));    
        for dd = 1:20
            margMat = zeros(1000, 2);
            for ss = 1:1000
                margMat(ss, :) = history{1, ss}(dd, :);
            end
            csvwrite(strcat('/Shared/ssrivastva/wasp/parafac/result/comp/sub5/res_cv_', num2str(cc), ...
                            '_sub_', num2str(kk), '_dim_', num2str(dd), '_k5.csv'), margMat);    
        end   
    end    
end

clear;

for cc = 1:10
    for kk=1:10
        load(strcat('/Shared/ssrivastva/wasp/parafac/result/comp/sub10/samp/res_cv_', ...
                    num2str(cc), '_sub_', num2str(kk), '_k10.mat'));    
        for dd = 1:20
            margMat = zeros(1000, 2);
            for ss = 1:1000
                margMat(ss, :) = history{1, ss}(dd, :);
            end
            csvwrite(strcat('/Shared/ssrivastva/wasp/parafac/result/comp/sub10/res_cv_', num2str(cc), ...
                            '_sub_', num2str(kk), '_dim_', num2str(dd), '_k10.csv'), margMat);    
        end   
    end    
end


