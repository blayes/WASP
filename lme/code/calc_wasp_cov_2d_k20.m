% dirp = '/Shared/ssrivastva/wasp/lme/result/wasp/'
function calc_wasp_cov_2d_k20(dd, nsub, ndim, dirp)

addpath('/opt/gurobi/6.5.1/linux64/matlab/');

rtime = zeros(ndim, 2);
grdsize = 60;
% calculate the pair-wise sq. euclidean distance between the atoms of subset
% posteriors and WASP atoms
for pp = 1:2 
    for dims = 1:ndim
        covs = {};
        for jj = 1:nsub
            covs{jj} = csvread(strcat(dirp, 'samp/joint/cov_cv_', num2str(dd), '_p_', ...
                                 num2str(pp),'_nsub_', num2str(jj), '_d_', ...
                                 num2str(dims), '_k20.csv'));         
        end
                
        subsetPost = {};    
        for jj = 1:nsub
            subsetPost{jj} = covs{jj}(randi([1 1000], 200, 1), :);
        end         

        lbd1 = min(cellfun(@(x) x(1), cellfun(@(x) min(x), subsetPost, ...
                                              'UniformOutput', false)));
        lbd2 = min(cellfun(@(x) x(2), cellfun(@(x) min(x), subsetPost, ...
                                              'UniformOutput', false)));   
        ubd1 = max(cellfun(@(x) x(1), cellfun(@(x) max(x), subsetPost, ...
                                              'UniformOutput', false)));
        ubd2 = max(cellfun(@(x) x(2), cellfun(@(x) max(x), subsetPost, ...
                                              'UniformOutput', false)));   

        [opostx, oposty] = meshgrid(linspace(lbd1, ubd1, grdsize), linspace(lbd2, ubd2, grdsize));
        overallPost = [opostx(:) oposty(:)]; % 
        
        distMatCell = {};
        
        m00 = diag(overallPost * overallPost');
        for ii = 1:nsub
            mm = diag(subsetPost{ii} * subsetPost{ii}');    
            mm1 = overallPost * subsetPost{ii}'; 
            distMatCell{ii} = bsxfun(@plus, bsxfun(@plus, -2 * mm1, mm'), m00);    
        end
        
        % constants
        K  = nsub;
        Ni = cell2mat(cellfun(@(x) size(x, 2), distMatCell, 'UniformOutput', false));
        N  = size(overallPost, 1);
        nx = N * (N+1);
        mx = K * N + N + 1;
        In = eye(N);
        En = ones(1, N);

        % Generate matrix A0.
        A0  = sparse([]);
        for p = 1:K
            cc = (1:N)';                  % terribly fast version of 
            idx = cc(:, ones(Ni(p), 1));  % repmat(In, 1, Ni(p)) / Ni(p)
            Rp  = In(:, idx(:)) / Ni(p);  % in 3 steps
            A0  = blkdiag(A0, Rp); 
        end
        cc = (1:N)';                  % terribly fast version of 
        idx = cc(:, ones(K, 1));      % repmat(-In, K, 1) 
        A00  = -In(idx(:), :);        % in 3 steps
        
        A0 = sparse([A00, A0]);
        b0 = zeros(size(A0, 1), 1);
        disp('done generating A ...');        
        
        % Generate matrix B from simplex constraints.
        B = sparse([]);
        for p = 0:(sum(Ni))
            B = blkdiag(B, En);
        end
        disp('done generating B ...');        
        
        % The hold matrix C.
        A = sparse([A0; B]);

        % Generate the right hand size vector b.
        b = sparse([zeros(K * N, 1); ones(sum(Ni) + 1, 1)]);
        
        % Generate the cost vector
        costCell = cellfun(@(x) x(:) / size(x, 2), distMatCell, 'UniformOutput', false);
        costVec = [zeros(size(overallPost, 1), 1); cell2mat(costCell(:))];
        
        c = sparse(costVec);
        
        tic;
        lpsol = callLpSolver('gurobi', A, b, c, 10000, 1e-10);
        rtime(dims, pp) = toc;
        
        [tmats, avec] = recoverSolution(lpsol, K, N, Ni);

        summ = [overallPost avec];
        csvwrite(strcat(dirp, 'joint/wasp_cov_cv_',  num2str(dd), '_p_', num2str(pp), '_d_', num2str(dims), '_k20.csv'), summ);
        
        disp(['done with sim ' num2str(dd) '...' ' p ' num2str(pp) '...' ' dim ' num2str(dims) '... ']);    
    end    
end

csvwrite(strcat(dirp, 'joint/cov_2d_times_cv_', num2str(dd),  '_k20.csv'), rtime);    

quit

