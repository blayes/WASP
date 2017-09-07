% FUNCTION:  [xsol, output] = callLpSolver(solver, Amat, bvec, ...
%                                          cvec, maxiters, tolx)
% PURPOSE: Call the Lp solver to solve the LP problem.
% 
%
function [xsol, output] = callLpSolver(solver, Amat, bvec, cvec, maxiters, tolx)

% Check the inputs.
if nargin < 3, error('At least three inputs are required!'); end
if nargin < 5, tolx = 1e-4; end
if nargin < 4, maxiters = 1000; end
if isempty(tolx), tolx = 1e-4; end
if isempty(maxiters), maxiters = 1000; end
if size(Amat, 1) ~= length(bvec), error('Inputs are inconsistent!'); end
if size(Amat, 2) ~= length(cvec), error('Inputs are inconsistent!'); end
xsol = []; output = [];
nx   = length(cvec);

%% Call the SeDuMi solver.
if strcmpi(solver, 'sedumi')
    pars.maxiter = maxiters;
    pars.eps     = tolx;
    Cones.l      = length(cvec);
    time3        = tic;
    [x_sedumi, y_sedumi, info] = sedumi(Amat, bvec, cvec, Cones, pars);
    x_sedumi = full(x_sedumi);
    xsol         = x_sedumi;
    output.time3        = toc(time3);
    output.info  = info;
    output.dual_sol = y_sedumi;
end

%% Call the Matlab LINPROG solver.
if strcmpi(solver, 'linprog')
    time5       = tic;
    opts        = optimset('Algorithm', 'interior-point', ...
                           'Display', 'iter', ...
                           'MaxIter', maxiters, 'TolX', tolx, 'TolFun', tolx);
    [xsol, fx2] = linprog(cvec, [], [], Amat, bvec, zeros(nx, 1), [], [], opts);
    output.time = toc(time5);
    output.fx   = fx2;
end

%% Call the SDPT3 solver.
if strcmpi(solver, 'sdpt3')
    Asdpt3 = spconvert(Amat);
    blk{1, 1} = 'l';
    blk{1, 2} = ones(1, nx);
    sdpt3opt  = struct('gaptol', tolx, 'maxit', maxiters);
    time4     = tic;
    %[X0, y0, Z0] = infeaspt(blk, Asdpt3, cvec, bvec); 
    %[fx_sdpt3, x_sdpt3, y_sdpt3, Z_sdpt3] = sqlp(blk, Asdpt3, cvec, bvec, sdpt3opt, X0, y0, Z0);
    [fx_sdpt3, x_sdpt3, y_sdpt3, Z_sdpt3] = sqlp(blk, Asdpt3, cvec, bvec, sdpt3opt);
    output.time4     = toc(time4);
    x_sdpt3   = x_sdpt3{:};
    xsol      = x_sdpt3;
    output.dual_sol = y_sdpt3;
    output.slacks   = Z_sdpt3;
    output.fx       = fx_sdpt3;
end    

%% Call Our Decopt solver.
if strcmpi(solver, 'decopt')
    
    % Set the parameters.
    param.MaxIters      = maxiters;
    param.Verbosity     = 2;
    param.RelTolX       = tolx;
    param.saveHistMode  = 0;
    param.Algorithm     = 3;
    param.InnerMaxIters = 20;
    param.adaptStepSize = 0;

    % Call the solver.
    proxLpPos    = @(x, gamma)( min( max(0, x - gamma*cvec), 1.0) );

    % User-define proximal-functions.
    proxOpers{1} = @(x, gamma, varargin)(proxLpPos(x, gamma));
    proxOpers{2} = @(x, gamma, varargin)(projL2norm(x, 1e-12));

    proxOpers{3} = @(x, varargin)( cvec'*x );
    proxOpers{4} = @(x, varargin)(0);
    
    % Generate an initial point.
    x0       = zeros(nx, 1);
    
    %% Call the solver with user-define prox-functions.
    time1 = tic;
    [xsol, out] = decoptSolver('UserDef', Amat, bvec, param, 'x0', x0, 'Prox', proxOpers, 'GammaFactor', 1.1);
    output.time = toc(time1);
    output.info = out;

end

%% Call the Gurobi solver.
if strcmpi(solver, 'gurobi')
    
    % Generate the LP model.
    time_g            = tic;
    model.A           = Amat;
    model.obj         = full(cvec);
    model.rhs         = full(bvec);
    model.modelsense  = 'min';
    model.sense       = '=';
    
    % Define the parameters.
    param.method      = 2;
    param.Presolve    = 2;
    param.Crossover   = 0;
    param.outputflag  = 1;
    
    % Call the solver.
    result            = gurobi(model, param);
    
    % Obtain the final results.
    output.result     = result;
    output.time       = toc(time_g);
    xsol              = result.x;
    
end

