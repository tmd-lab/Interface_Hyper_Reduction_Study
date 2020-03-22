function [R, dRdX, dRda, varargout] = WJRES_STR_WTH(X, pars, K, Fd, Fs, L, Txyntxyn, Qxyntxyn, Fxynfun, Nps, varargin)
%WJRESFUN Residual Function for Whole Joint Models
% USAGE:
%  [R, dRdX, dRda] = WJRESFUN(X, K, Fd, Fs, L, Fxfun, Fyfun, Fnfun);
% INPUTS:
%   X       :
%   K       :
% OUTPUTS:
%   R
%   dRdX
%   dRda
    if nargin>=13
        RCvec = varargin{1};
        RX = RCvec.*X;
    else
        RCvec = ones(size(X));
        RX = X;
    end
    
    alpha = RX(end);
    uxyntxyn = [Qxyntxyn(1:6:end,:)*RX(1:end-1), Qxyntxyn(2:6:end,:)*RX(1:end-1), Qxyntxyn(3:6:end,:)*RX(1:end-1), ...
        Qxyntxyn(4:6:end,:)*RX(1:end-1) Qxyntxyn(5:6:end,:)*RX(1:end-1) Qxyntxyn(6:6:end,:)*RX(1:end-1)];
    if nargout<4
        [Fxyntxyn, dFxyntxyn] = Fxynfun(uxyntxyn, pars);
    else        
        [Fxyntxyn, dFxyntxyn, dFxyndpar] = Fxynfun(uxyntxyn, pars);
    end
	if size(dFxyntxyn,3)==1
       tmp = dFxyntxyn;
       dFxyntxyn = zeros(size(dFxyntxyn,1),size(dFxyntxyn,2),3);
       dFxyntxyn(:,1,1) = tmp(:,1);
       dFxyntxyn(:,2,2) = tmp(:,2);
       dFxyntxyn(:,3,3) = tmp(:,3); % Np x uxyn x fxyn
    end
	dFxyntxyn = num2cell(permute(dFxyntxyn,[3,2,1]), [2,1]);
    
    Ng = size(L,1) - Nps*6;
    Fnl = Txyntxyn*reshape(Fxyntxyn', Nps*6, 1);
%     dFnl = Txyn*diag(reshape(dFxyn', Nps*3, 1))*Qxyn;
	dFnl = Txyntxyn*blkdiag(dFxyntxyn{:})*Qxyntxyn;
    
    R = K*RX(1:end-1) + Fnl - Fs - alpha*Fd;
    dRdX = K.*RCvec(1:end-1) + dFnl.*RCvec(1:end-1);
    dRda = -Fd*RCvec(end);
    
	if nargin>=14
        R = varargin{2}.*R;
        dRdX = varargin{2}.*dRdX;
        dRda = varargin{2}.*dRda;
    end
    if nargout>=4
        varargout{1} = Txyntxyn*reshape(permute(dFxyndpar, [3, 1, 2]), Nps*6, length(pars));
    end
end
