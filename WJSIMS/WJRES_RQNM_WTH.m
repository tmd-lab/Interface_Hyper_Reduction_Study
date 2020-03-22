function [R, dRdX, dRdq, varargout] = WJRES_RQNM_WTH(X, pars, prev, K, M, Xstat, dXsdp, Fs, L, Txyntxyn, Qxyntxyn, Fxynfun, Nps, varargin)
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
    if nargin>=14
        RCvec = varargin{1};
        RX = RCvec.*X;
    else
        RCvec = ones(size(X));
        RX = X;
    end
    
    q = RX(end);
    lam = RX(end-1);
    
    uxyntxyn = [Qxyntxyn(1:6:end,:)*RX(1:end-2), Qxyntxyn(2:6:end,:)*RX(1:end-2), Qxyntxyn(3:6:end,:)*RX(1:end-2), ...
        Qxyntxyn(4:6:end,:)*RX(1:end-2) Qxyntxyn(5:6:end,:)*RX(1:end-2) Qxyntxyn(6:6:end,:)*RX(1:end-2)];
    if nargout<4
        [Fxyntxyn, dFxyntxyn] = Fxynfun(uxyntxyn-prev.uxyntxyn, pars);
        Fxyntxyn = Fxyntxyn + prev.Fxyntxyn;
    else
        [Fxyntxyn, dFxyntxyn, dFxyntxyndpar] = Fxynfun(uxyntxyn-prev.uxyntxyn, pars);
        Fxyntxyn = Fxyntxyn + prev.Fxyntxyn;
        dFxyntxyndpar = dFxyntxyndpar + prev.dFxyntxyndpar;
    end
    if size(dFxyntxyn,3)==1
       tmp = dFxyntxyn;
       dFxyntxyn = zeros(size(dFxyntxyn,1),size(dFxyntxyn,2),3);
       dFxyntxyn(:,1,1) = tmp(:,1);
       dFxyntxyn(:,2,2) = tmp(:,2);
       dFxyntxyn(:,3,3) = tmp(:,3); % Np x uxyn x fxyn
    end
	dFxyntxyn = num2cell(permute(dFxyntxyn,[2,3,1]), [2,1]);
    
    Ng = size(L,1) - Nps*6;
    Fnl = Txyntxyn*reshape(Fxyntxyn', Nps*6, 1);
%     dFnl = Txyn*diag(reshape(dFxyn', Nps*3, 1))*Qxyn;
    dFnl = Txyntxyn*blkdiag(dFxyntxyn{:})*Qxyntxyn;
    
    R = [K*RX(1:end-2) + Fnl - Fs - lam*M*(RX(1:end-2)-Xstat);
        (RX(1:end-2)-Xstat)'*M*(RX(1:end-2)-Xstat)-q^2];
    dRdX = [(K - lam*M + dFnl).*RCvec(1:end-2), -M*(RX(1:end-2)-Xstat);
        2*(RX(1:end-2)-Xstat)'*M, 0];
    dRdq = [zeros(length(RX)-2,1); -2*q];
    
%     R = [K*RX(1:end-2) + Fnl - Fs - lam*M*(RX(1:end-2)-Xstat);
%         (RX(1:end-2)-Xstat)'*(K*RX(1:end-2) + Fnl - Fs) - lam*q^2];
%     dRdX = [(K - lam*M + dFnl).*RCvec(1:end-2), -M*(RX(1:end-2)-Xstat);
%         RX(1:end-2)'*K+Fnl'-Fs'+(RX(1:end-2)-Xstat)'*(K+dFnl), -q^2];
%     dRdq = [zeros(length(RX)-2,1); -2*lam*q];
    
	if nargin>=15
        R = varargin{2}.*R;
        dRdX = varargin{2}.*dRdX;
        dRdq = varargin{2}.*dRdq;
    end
    if nargout>=4
        varargout{1} = zeros(Nps*6, length(pars));
        for i=1:length(pars)
            varargout{1}(:, i) = reshape(squeeze(dFxyntxyndpar(:,i,:))', Nps*6, 1);
        end
%         varargout{1} = reshape(permute(dFxyndpar, [3 1 2]), 3*Nps, length(pars))-(reshape(dFxyn', Nps*3, 1)).*Qxyn*dXsdp+prev.dFxyndXdpar;
        varargout{1} = reshape(permute(dFxyntxyndpar, [3 1 2]), Nps*6, length(pars))-blkdiag(dFxyntxyn{:})*Qxyntxyn*dXsdp+prev.dFxyntxyndXdpar;
        varargout{1} = [Txyntxyn*varargout{1}+lam*M*dXsdp; -2*(RX(1:end-2)-Xstat)'*M*dXsdp];
    end
end