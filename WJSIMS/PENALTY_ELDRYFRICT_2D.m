function [Fxyn, varargout] = PENALTY_ELDRYFRICT_2D(uxyn, pars, opt)
%JENKINS friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = PENALTY_ELDRYFRICT(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [mu; Kt/A; Kn/A]
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Fxyn - Force vectors, rows are partches, columns are x,y,z
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Penalty Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = max(Kn.*uxyn(:,3),0);
    coni = (Fxyn(:,3)~=0);
    
    % Parameters - X
    Mu = opt.x.T{1}*lpars;
    Fs = Mu.*Fxyn(:,3);
    Ktx = opt.x.T{2}*lpars;
    
	% Parameters - Y
    Kty = opt.y.T{2}*lpars;
    
    Fxsp = Ktx.*uxyn(:,1);
    Fysp = Kty.*uxyn(:,2);
    Ftsp = sqrt(Fxsp+Fysp);
    
    Fxn = Fxsp./sqrt(Fxsp.^2+Fysp.^2);  Fxn(~isfinite(Fxn)) = 0;
    Fyn = Fysp./sqrt(Fxsp.^2+Fysp.^2);  Fyn(~isfinite(Fyn)) = 0;
    
    Fxyn(:,1) = (Fxsp.*(Ftsp < Fs) + Fs.*Fxn.*(Ftsp >= Fs)).*coni;
    Fxyn(:,2) = (Fysp.*(Ftsp < Fs) + Fs.*Fyn.*(Ftsp >= Fs)).*coni;
    
    sux = sign(real(uxyn(:,1)));
    ux = sux.*uxyn(:,1);

    suy = sign(real(uxyn(:,2)));
    uy = suy.*uxyn(:,2);    
    %% Jacobians
    if nargout>=2
        
        %this is [partial F_x/partial ux, partial Fy/partial uy, partial
        %Fz/partial uz]
        varargout{1} = zeros(size(Fxyn,1),3,3);
        % X Force
        varargout{1}(:,1,1) = (Ktx.*(Ftsp < Fs) + (Fs.*Fyn.^2.*Ktx./Ftsp).*(Ftsp >= Fs)).*coni;
        varargout{1}(:,2,1) = -(Fs.*Fxn.*Fyn.*Kty./Ftsp).*(Ftsp >= Fs).*coni;
        varargout{1}(:,3,1) = Mu.*Kn.*Fxn.*(Ftsp >= Fs).*coni;
        % Y Force
        varargout{1}(:,1,2) = -(Fs.*Fxn.*Fyn.*Ktx./Ftsp).*(Ftsp >= Fs).*coni;
        varargout{1}(:,2,2) = (Kty.*(Ftsp < Fs) + (Fs.*Fxn.^2.*Kty./Ftsp).*(Ftsp >= Fs)).*coni;
        varargout{1}(:,3,2) = Mu.*Kn.*Fyn.*(Ftsp >= Fs).*coni;
        % N Force
        varargout{1}(:,3,3) = Kn.*coni;
        
%         varargout{1} = [Kt.*(Kt.*ux < Fs), kt.*(kt.*uy < fs), Kn].*coni;
    end
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        %Jenkins - X
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*Fxyn(:,3).*sux.*(Ftsp >= Fs), ...
            opt.x.T{2}(:, 2).*sux.*ux.*(Ftsp < Fs),...
            2*opt.x.T{1}(:,1).*((opt.n.T{1}(:,3)*lpars(3)).*uxyn(:,3)).*sux.*(Ftsp >= Fs)];
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*coni.*dlparsdpars';  % Check        
        
        %Jenkins - Y
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*suy.*(Ftsp >= Fs), ...
            opt.y.T{2}(:, 2).*suy.*uy.*(Ftsp < Fs),...
            2*opt.y.T{1}(:,1).*((opt.n.T{1}(:,3)*lpars(3)).*uxyn(:,3)).*suy.*(Ftsp >= Fs)];
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*coni.*dlparsdpars';  % Check        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyn, 1), length(pars)-1) (uxyn(:,3).*opt.n.T{1}(:, 3)).*coni].*dlparsdpars';
    end
end