function [Fxyntxyn, varargout] = PENALTY_ELDRYFRICT_2D_WTH(uxyntxyn, pars, opt)
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
%   Fxyn - Force vectors, rows are partches, columns are x,y,z,tx,ty,tz
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    %% Forces    
    Fxyntxyn = zeros(size(uxyntxyn));
    
    % Penalty Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyntxyn(:, 3) = max(Kn.*uxyntxyn(:,3),0);
    coni = (Fxyntxyn(:,3)~=0);
    
    % Parameters - X
    Mu = opt.x.T{1}*lpars;
    Fs = Mu.*Fxyntxyn(:,3);
    Ktx = opt.x.T{2}*lpars;
    
	% Parameters - Y
    Kty = opt.y.T{2}*lpars;
    
    % Planar Forces
    Fxsp = Ktx.*uxyntxyn(:,1);
    Fysp = Kty.*uxyntxyn(:,2);
    Ftsp = sqrt(Fxsp.^2+Fysp.^2);
    
    Fxn = Fxsp./Ftsp;  Fxn(~isfinite(Fxn)) = 0;
    Fyn = Fysp./Ftsp;  Fyn(~isfinite(Fyn)) = 0;
    
    Fxyntxyn(:,1) = (Fxsp.*(Ftsp < Fs) + Fs.*Fxn.*(Ftsp >= Fs)).*coni;
    Fxyntxyn(:,2) = (Fysp.*(Ftsp < Fs) + Fs.*Fyn.*(Ftsp >= Fs)).*coni;
    
    sux = sign(real(uxyntxyn(:,1)));
    ux = sux.*uxyntxyn(:,1);

    suy = sign(real(uxyntxyn(:,2)));
    uy = suy.*uxyntxyn(:,2);    
    
    % Rotational moments
%     Ktxx = opt.tx.T{1}*lpars; Ktxy = opt.tx.T{2}*lpars;
%     Ktyx = opt.ty.T{1}*lpars; Ktyy = opt.ty.T{2}*lpars;    
%     Ktnn = opt.tn.T{1}*lpars;
%     Fxyn(:,4) = (Ktxx.*uxyntxyn(:,4) + Ktxy.*uxyntxyn(:,5)).*coni;
%     Fxyn(:,5) = (Ktyx.*uxyntxyn(:,4) + Ktyy.*uxyntxyn(:,5)).*coni;
%     Fxyn(:,6) = Ktnn.*uxyntxyn(:,6).*coni;
%     
    Ixx = opt.tx.T{1}*ones(size(lpars)); Ixy = opt.tx.T{2}*ones(size(lpars));
    Iyx = opt.ty.T{1}*ones(size(lpars)); Iyy = opt.ty.T{2}*ones(size(lpars));
    Izz = opt.tn.T{1}*ones(size(lpars));
    Fxyntxyn(:,4) = (Ixx.*lpars(2).*uxyntxyn(:,4)+Ixy.*lpars(2).*uxyntxyn(:,5)).*(Ftsp<Fs).*coni;
    Fxyntxyn(:,5) = (Iyx.*lpars(2).*uxyntxyn(:,4)+Iyy.*lpars(2).*uxyntxyn(:,5)).*(Ftsp<Fs).*coni;
    Fxyntxyn(:,6) = Izz.*lpars(3).*uxyntxyn(:,6).*coni;    %% Jacobians
    if nargout>=2
        
        %this is [partial F_x/partial ux, partial Fy/partial uy, partial
        %Fz/partial uz]
        varargout{1} = zeros(size(Fxyntxyn,1),6,6);
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
        % TX Moment
        varargout{1}(:,4,4) = Ixx*lpars(2).*(Ftsp<Fs).*coni;
        varargout{1}(:,5,4) = Ixy*lpars(2).*(Ftsp<Fs).*coni;
        % TY Moment
        varargout{1}(:,4,5) = Iyx*lpars(2).*(Ftsp<Fs).*coni;
        varargout{1}(:,5,5) = Iyy*lpars(2).*(Ftsp<Fs).*coni;
        % TN Moment
        varargout{1}(:,6,6) = Izz*lpars(3).*(Ftsp<Fs).*coni;
        
%         varargout{1} = [Kt.*(Kt.*ux < Fs), kt.*(kt.*uy < fs), Kn].*coni;
    end
    if nargout>=3
        varargout{2} = zeros(size(Fxyntxyn,1), length(lpars), 6);
        
        %Jenkins - X
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*Fxyntxyn(:,3).*sux.*(Ftsp >= Fs), ...
            opt.x.T{2}(:, 2).*sux.*ux.*(Ftsp < Fs),...
            2*opt.x.T{1}(:,1).*((opt.n.T{1}(:,3)*lpars(3)).*uxyntxyn(:,3)).*sux.*(Ftsp >= Fs)];
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*coni.*dlparsdpars';  % Check        
        
        %Jenkins - Y
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*suy.*(Ftsp >= Fs), ...
            opt.y.T{2}(:, 2).*suy.*uy.*(Ftsp < Fs),...
            2*opt.y.T{1}(:,1).*((opt.n.T{1}(:,3)*lpars(3)).*uxyntxyn(:,3)).*suy.*(Ftsp >= Fs)];
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*coni.*dlparsdpars';  % Check        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyntxyn, 1), length(pars)-1) (uxyntxyn(:,3).*opt.n.T{1}(:, 3)).*coni].*dlparsdpars';
    end
end