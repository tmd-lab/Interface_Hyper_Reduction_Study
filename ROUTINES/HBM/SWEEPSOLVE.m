function [X] = SWEEPSOLVE(func, X0, W, ematfun, opt)
%SWEEPSOLVE Sweep system between two frequencies using tangent predictor
% USAGE:
%   [X] = SWEEPSOLVE(func, X0, W, ematfun, opt)
% INPUTS:
%   func        : Function handle
%   X0          : 
%   W           :      
%   ematfun     : Function handle
%   opt         : Optimoptions
% OUTPUTS:
%   X           :
    % Initial Correction
    [X0, ~, efl, ~, Jf] = fsolve(@(X) func([X; W(1)]), X0, opt);
    if efl <= 0
        error('Initial Correction Unsuccessful - Quitting!');
    end
    disp('-------------------------------------------------');
    disp('INITIAL SOLUTION CONVERGED');
    disp('-------------------------------------------------');
    [~, dEdw] = ematfun(W(1));
    dRdw = -Jf\dEdw*X0;
    
    % Sweep Solve
    Nw = length(W);
    X = zeros(length(X0)+1, Nw);
    X(:, 1) = [X0; W(1)];
    for i=2:Nw
        X0 = X(1:end-1,i-1) + dRdw*(W(i)-W(i-1));  % First order predictor
        [X0, ~, efl, ~, Jf] = fsolve(@(X) func([X; W(i)]), X0, opt);
        if efl <= 0
            disp(['Unsuccessful Step ' num2str(i)]);
            break
        end
        X(:, i) = [X0; W(i)];
        [~, dEdw] = ematfun(W(i));
        dRdw = -Jf\dEdw*X0;
        disp('-------------------------------------------------');
        fprintf('DONE %d/%d\n',i,Nw);
        disp('-------------------------------------------------');
    end
end

