function [Ptx,Pty,Pn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun] ...
        = JENKINS_P(uxyn,pn0,pars,uxyntxynp) 
%JENKINS_P 3D Jenkins model Pressure and Jacobians
% USAGE:
%	[Ptx,Pty,Pn,Jttx,Jtnx,Jtty,Jtny,Jnn] = JENKINS_P(uxyn,pn0,pars,uxyntxynp);
% INPUTS:
%   uxyn	: (Npx3) 2 tangential and 1 normal displacement vectors 
%   pn0         : Static Normal pressure
%   pars        : Structure containing:
%       mu	: coefficient of friction (scalar)
%       ktx	: x tangential stiffness (scalar)
%       kty	: y tangential stiffness (scalar)
%       kn	: normal stiffness (scalar)
%   uxyntxynp	: (Npx6) vectors of displacements & tractions at
%   		  start state
% OUTPUTS:
%   Ptx         : Npx1 x tangential traction
%   Pty         : Npx1 y tangential traction
%   Pn          : Npx1 normal pressure
%   dtxdux      : Npx1 dtx/dux
%   dtxduy      : Npx1 dtx/duy
%   dtxdun      : Npx1 dtx/dun    
%   dtydux      : Npx1 dty/dux
%   dtyduy      : Npx1 dty/duy
%   dtydun      : Npx1 dty/dun        
%   dtndun      : Npx1 dtn/dun

    mu 	= pars.mu;
    ktx = pars.ktx;
    kty = pars.kty;
    kn  = pars.kn;
    
    ux  = uxyn(:,1);
    uy  = uxyn(:,2);
    un  = uxyn(:,3);
    uxp = uxyntxynp(:,1);
    uyp = uxyntxynp(:,2);
    unp = uxyntxynp(:,3);
    txp = uxyntxynp(:,4);
    typ = uxyntxynp(:,5);
    tnp = uxyntxynp(:,6);
    
    Np = size(un,1);    
    
    prevsep    = find((pn0+tnp)<=0);
    prevstick  = find(sqrt(txp.^2+typ.^2)<=abs(mu*(pn0+tnp)));
    prevslip   = setdiff(setdiff(1:Np,prevsep),prevstick);

    pn = max(pn0 + kn*un,0);    % Normal force with penalty spring        
    pxstick = ktx.*(ux-uxp);          % Stick traction prediction x    
%     pxstick(prevstick) = pxstick(prevstick) + txp(prevstick); % Memory preserved if it was sticking at turning point
    pxstick = pxstick + txp;
    pystick = kty.*(uy-uyp);          % Stick traction prediction y
%     pystick(prevstick) = pystick(prevstick) + typ(prevstick); % Memory preserved if it was sticking at turning point
    pystick = pystick + typ;
    pslip   = mu.*pn;                             % Slip traction
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uxynorm = sqrt((pxstick).^2+(pystick).^2);
    uxynorm(abs(uxynorm)<eps) = 1.0;
    ux_normed = (pxstick)./uxynorm;  ux_normed(abs(ux_normed)<eps) = 0.0;
    uy_normed = (pystick)./uxynorm;  uy_normed(abs(uy_normed)<eps) = 0.0;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sep = find(pn==0); 		% Separated cases
    cnt = setdiff(1:Np, sep); 	% Contacted cases    
    
    slipped = find(sqrt(pxstick.^2+pystick.^2)>abs(pslip));
    slipped = setdiff(slipped,sep);
    stuck = setdiff(setdiff(1:Np,slipped),sep);
    
    Ptx    	= zeros(Np,1);
    Pty    	= zeros(Np,1);
    Pn     	= zeros(Np,1);
    
    dtxdux 	= zeros(Np,1);	dtxduxstick = zeros(Np,1);	
    dtxduy 	= zeros(Np,1);
    dtxdun 	= zeros(Np,1);
    
    dtydux 	= zeros(Np,1);
    dtyduy 	= zeros(Np,1);	dtyduystick = zeros(Np,1);
    dtydun 	= zeros(Np,1);
    
    dtndun 	= zeros(Np,1);
    
    % Normal Pressure & Jacobian
    Pn = pn - pn0;
    dtndun(cnt) = kn;
    
    % X direction traction & Jacobians
    pxslip = pslip.*ux_normed; 	% Slip traction x %%%%%%%%%%%%%%%%
    dtxduxstick(:) 	= ktx;
    dtyduystick(:) 	= kty;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dtxduxslip      = mu.*pn.*(pystick).^2.*dtxduxstick./(uxynorm).^3.0;
    dtxduyslip      = -mu.*pn.*(pxstick.*pystick).*dtyduystick./(uxynorm).^3.0;
    dtxdunslip      = mu.*kn.*ux_normed;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ptx(stuck) 		= pxstick(stuck);
    Ptx(slipped) 	= pxslip(slipped);
    Ptx(sep) 		= 0.0;

    dtxdux(stuck)	= dtxduxstick(stuck);
    dtxdux(slipped)	= dtxduxslip(slipped);
    dtxdux(sep)		= 0.0;
    
    dtxduy(stuck)	= 0.0;
    dtxduy(slipped)	= dtxduyslip(slipped);
    dtxduy(sep)		= 0.0;

    dtxdun(stuck)	= 0.0;
    dtxdun(slipped)	= dtxdunslip(slipped);
    dtxdun(sep)		= 0.0;
    
    % Y direction traction & Jacobian
    pyslip = pslip.*uy_normed; 	% Slip traction y %%%%%%%%%%%%%%%%
    dtyduystick(:) 	= kty;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dtyduxslip      = -mu.*pn.*(pxstick.*pystick).*dtxduxstick./(uxynorm).^3.0;
    dtyduyslip      = mu.*pn.*(pxstick).^2.*dtyduystick./(uxynorm).^3.0;
    dtydunslip      = mu.*kn.*uy_normed;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pty(stuck) 		= pystick(stuck);
    Pty(slipped) 	= pyslip(slipped);
    Pty(sep) 		= 0.0;
    
    dtydux(stuck)	= 0.0;
    dtydux(slipped)	= dtyduxslip(slipped);
    dtydux(sep)		= 0.0;
    
    dtyduy(stuck)	= dtyduystick(stuck);
    dtyduy(slipped)	= dtyduyslip(slipped);
    dtyduy(sep)		= 0.0;
    
    dtydun(stuck)	= 0.0;
    dtydun(slipped)	= dtydunslip(slipped);
    dtydun(sep)		= 0.0;

%     disp([num2str(length(sep)/Np*100) ' percentage separated'])
%     disp([num2str(length(slipped)/Np*100) ' percentage slipped'])
%     disp([num2str(length(stuck)/Np*100) ' percentage stuck'])
end