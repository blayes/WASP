% FUNCTION: [tOptSol, aOptSol] = recoverSolution(xOptSol, nsubs, N, Ni)
% PURPOSE:  Recover the original solution.
%
function [tOptSol, aOptSol] = recoverSolution(xOptSol, nsubs, N, Ni)

    % Define this soft-thresholding operator to remove small elements.
    softThresOper = @(x, t)(sign(x).*max(abs(x) - t, 0));

    % Recover the solution.
    aOptSol = xOptSol(1:N, 1);
    xRest   = xOptSol(N+1:end);
    for p = 1:nsubs
        tOptSol{p} = reshape( xRest(1:N*Ni(p)), N, Ni(p))/Ni(p);
        xRest      = xRest(N*Ni(p)+1:end);
        tOptSol{p} = softThresOper(tOptSol{p}, 1e-10);
    end

end
% @END ...
