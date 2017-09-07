%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Replicates simulation study of Dunson and Xing (2009) 
% 
% http://www.tandfonline.com/doi/abs/10.1198/jasa.2009.tm08439#.Uxpc6Nww_0A
%
% based on a version by Jing Zhou of UNC, Biostatistics
% modified by SS 05/10/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
rng(12345)

% -- global parameters -- %

% N = #sample
% q = #dimensions
% d = #categories in each dim.
% rep = #replications
N = 100000; q = 20; d = 2; rep = 10;

F1 = 0.5 * ones(q, d);    % category probabilities
F2 = F1; % category probabilities of response in sub-populations

F1(2,:) = [0.20 0.80]; F1(4,:) = [0.25 0.75]; F1(12,:) = [0.80 0.20]; F1(14,:) = [0.75 0.25];
F2(2,:) = [0.80 0.20]; F2(4,:) = [0.75 0.25]; F2(12,:) = [0.20 0.80]; F2(14,:) = [0.25 0.75];   

Yn_big = zeros(N, q, rep);
Si = ones(N, 1);

locs1 = randsample(1:N/2, 0.2*(N/2)); locs2 = randsample(N/2+1:N, 0.2*(N/2));
locs = [locs1 locs2];
Si(locs) = 0;

for g = 1:rep
    for j = 1:q
        Yn_big(Si==1,j,g) = mnrnd(1,F1(j,:),sum(Si==1))*(1:d)';
        Yn_big(Si==0,j,g) = mnrnd(1,F2(j,:),sum(Si==0))*(1:d)';
    end
end

save('/Shared/ssrivastva/wasp/parafac/data/parafac_full_data.mat', 'Yn_big', 'Si', 'locs', 'locs1', 'locs2', 'F1', 'F2', 'N', 'q', 'd', 'rep');    

clear;clc;
rng(12345)

load('/Shared/ssrivastva/wasp/parafac/data/parafac_full_data.mat');
% -- global parameters -- %

nsub = 5;

si1 = find(Si == 1);
si0 = find(Si == 0);

for cc = 1:rep
    si1Part = randi(nsub, length(si1), 1);
    si0Part = randi(nsub, length(si0), 1);    
    for kk = 1:nsub
        n1s = si1(si1Part == kk);
        n0s = si0(si0Part == kk);        
        partMat = Yn_big([n1s; n0s], :, cc);
        disp(['cv: ' num2str(cc) ' dim: ' num2str(size(partMat))]);
        save(strcat('/Shared/ssrivastva/wasp/parafac/data/sub5/parafac_cv_', ...
                    num2str(cc), '_sub_', num2str(kk), '_k5.mat'), 'partMat', 'n1s', 'n0s', 'si1Part', 'si0Part');            
    end    
end

clear;clc;
rng(12345)

load('/Shared/ssrivastva/wasp/parafac/data/parafac_full_data.mat');

% -- global parameters -- %

nsub = 10;

si1 = find(Si == 1);
si0 = find(Si == 0);

for cc = 1:rep
    si1Part = randi(nsub, length(si1), 1);
    si0Part = randi(nsub, length(si0), 1);    
    for kk = 1:nsub
        n1s = si1(si1Part == kk);
        n0s = si0(si0Part == kk);        
        partMat = Yn_big([n1s; n0s], :, cc);
        disp(['cv: ' num2str(cc) ' dim: ' num2str(size(partMat))]);        
        save(strcat('/Shared/ssrivastva/wasp/parafac/data/sub10/parafac_cv_', ...
                    num2str(cc), '_sub_', num2str(kk), '_k10.mat'), 'partMat', 'n1s', 'n0s', 'si1Part', 'si0Part');            
    end    
end

% $$$ % --- test code --- %
% $$$ [history, tend] = parafac_dx_sub(Yn_big(1:100, :, 1), cats, 20, 100000, 2000, 1000, 5); 
% $$$ [fhistory, tend] = parafac_dx_com(Yn_big(500,:,1), cats, 20, 2000, 1000, 5);           
% $$$ 
% $$$ load full/res_2.mat
% $$$ 
% $$$ full=history;
% $$$ fmargMat = zeros(1000, 2);
% $$$ for ss = 1:1000
% $$$     fmargMat(ss, :) = full{1, ss}(4, :);
% $$$ end
% $$$ 
% $$$ load sub10/samp/res_cv_2_sub_2_k10_6_27.mat
% $$$ 
% $$$ part = history;
% $$$ margMat = zeros(1000, 2);
% $$$ for ss = 1:1000
% $$$     margMat(ss, :) = part{1, ss}(4, :);
% $$$ end
% $$$ 
% $$$ hold off;
% $$$ plot(margMat(:, 1), 'Color', [1 0 1] );
% $$$ hold on;
% $$$ plot(fmargMat(:, 1));
% $$$ % 
% $$$ 
% $$$ hold off;
% $$$ hist(fmargMat(:, 1))
% $$$ hold on;
% $$$ hist(margMat(:, 1))
% $$$ 
% $$$ subplot(2, 1, 1)
% $$$ 
% $$$ subplot(2, 1, 2)
% $$$ plot(margMat(:, 1), 'Color', [1 0 1])
