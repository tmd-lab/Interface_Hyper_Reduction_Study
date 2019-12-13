function [Ptx,Pty,Pn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun] ...
        = ROUGHCONT_P(uxyn,pn0,cpars,uxyntxynp) 
%JENKINS_P 3D Rough Contact model Pressure and Jacobians
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

    mu 	= cpars.mu;
    ktx = cpars.ktx;
    kty = cpars.kty;
    kn  = cpars.kn;
    
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
    Ptx    	= zeros(Np,1);
    Pty    	= zeros(Np,1);
    Pn     	= zeros(Np,1);
    
    dtxdux 	= zeros(Np,1);	dtxduxstick = zeros(Np,1);	
    dtxduy 	= zeros(Np,1);  dtxdun 	= zeros(Np,1);
    
    dtydux 	= zeros(Np,1);  dtyduystick = zeros(Np,1);
    dtyduy 	= zeros(Np,1);	dtydun 	= zeros(Np,1);
    
    dtndun 	= zeros(Np,1);
    pn      = zeros(Np,1);
    if length(pn0)==1
        pn0 = ones(Np,1)*pn0;
    end
    
    prevsep    = find((pn0+tnp)<=0);
    prevstick  = find(sqrt(txp.^2+typ.^2)<=abs(mu.*(pn0+tnp)));
    prevslip   = setdiff(setdiff(1:Np,prevsep),prevstick);

    % Normal force with penalty condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sep = find(un<=0);  cnt = setdiff((1:Np)', sep); 	% Separation & Contact
    pn(cnt) = pn0(cnt) + kn(cnt, 1).*un(cnt).^kn(cnt, 2);
    pslip   = mu.*pn;                               % Slip traction
    reg1 = intersect(cnt, find(un<=cpars.kthr));    % Operation regime 1
    reg2 = setdiff(cnt, find(un<=cpars.kthr));      % Operation regime 2
    unc = un;  unc(sep) = 0;  % un but zero where in separation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Stuck Predictors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pxstick = zeros(Np, 1);
    pxstick(reg1) = (ktx(reg1,1).*unc(reg1).^ktx(reg1,2)).*(ux(reg1)-uxp(reg1)) + ...
        (ktx(reg1,1).*ktx(reg1,2).*unc(reg1).^(ktx(reg1,2)-1)).*ux(reg1).*(un(reg1)-unp(reg1));
    pxstick(reg2) = (ktx(reg2,3).*unc(reg2).^ktx(reg2,4)).*(ux(reg2)-uxp(reg2)) + ...
        (ktx(reg2,3).*ktx(reg2,4).*unc(reg2).^(ktx(reg2,4)-1)).*ux(reg2).*(un(reg2)-unp(reg2));
%     pxstick(prevstick) = pxstick(prevstick) + txp(prevstick); % Memory preserved if it was sticking at turning point
    pxstick = pxstick + txp;
    
    dtxduxstick(reg1) 	= ktx(reg1,1).*un(reg1).^ktx(reg1,2).*(1.0+ktx(reg1,2).*(1.0-unp(reg1)./un(reg1)));
    dtxduxstick(reg2) 	= ktx(reg2,3).*un(reg2).^ktx(reg2,4).*(1.0+ktx(reg2,4).*(1.0-unp(reg2)./un(reg2)));
    
    dtxdunstick(reg1)   = ktx(reg1,1).*ktx(reg1,2).*unc(reg1).^(ktx(reg1,2)-1).*(ux(reg1).*(2+(ktx(reg1,2)-1).*(1-unp(reg1)./unc(reg1))-uxp(reg1)));
    dtxdunstick(reg2)   = ktx(reg2,3).*ktx(reg2,4).*unc(reg2).^(ktx(reg2,4)-1).*(ux(reg2).*(2+(ktx(reg2,4)-1).*(1-unp(reg2)./unc(reg2))-uxp(reg2)));
    
    pystick = zeros(Np, 1);
    pystick(reg1) = (kty(reg1,1).*unc(reg1).^kty(reg1,2)).*(uy(reg1)-uyp(reg1)) + ...
        (kty(reg1,1).*kty(reg1,2).*unc(reg1).^(kty(reg1,2)-1)).*uy(reg1).*(un(reg1)-unp(reg1));
    pystick(reg2) = (kty(reg2,3).*unc(reg2).^kty(reg2,4)).*(uy(reg2)-uyp(reg2)) + ...
        (kty(reg2,3).*kty(reg2,4).*unc(reg2).^(kty(reg2,4)-1)).*uy(reg2).*(un(reg2)-unp(reg2));
%     pystick(prevstick) = pystick(prevstick) + typ(prevstick); % Memory preserved if it was sticking at turning point
    pystick = pystick + typ;
    
    dtyduystick(reg1) 	= kty(reg1,1).*un(reg1).^kty(reg1,2).*(1.0+kty(reg1,2).*(1.0-unp(reg1)./un(reg1)));
    dtyduystick(reg2) 	= kty(reg2,3).*un(reg2).^kty(reg2,4).*(1.0+kty(reg2,4).*(1.0-unp(reg2)./un(reg2)));
    
    dtydunstick(reg1)   = kty(reg1,1).*kty(reg1,2).*unc(reg1).^(kty(reg1,2)-1).*(uy(reg1).*(2+(kty(reg1,2)-1).*(1-unp(reg1)./unc(reg1))-uyp(reg1)));
    dtydunstick(reg2)   = kty(reg2,3).*kty(reg2,4).*unc(reg2).^(kty(reg2,4)-1).*(uy(reg2).*(2+(kty(reg2,4)-1).*(1-unp(reg2)./unc(reg2))-uyp(reg2)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Setup Slip directions if need be %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uxynorm = sqrt(pxstick.^2+pystick.^2);
    uxynorm(abs(uxynorm)<eps) = 1.0;
    ux_normed = pxstick./uxynorm;  ux_normed(abs(ux_normed)<eps) = 0.0;
    uy_normed = pystick./uxynorm;  uy_normed(abs(uy_normed)<eps) = 0.0;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine Stuck & Slipped Cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slipped = find(sqrt(pxstick.^2+pystick.^2)>abs(pslip));
    slipped = setdiff(slipped,sep);
    stuck = setdiff(setdiff(1:Np,slipped),sep);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Normal Traction & Jacobian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pn = pn - pn0;
    dtndun(cnt) = kn(cnt, 1).*kn(cnt, 2).*unc(cnt).^(kn(cnt, 2)-1.0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % X-Tangential tractions & jacobians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pxslip = pslip.*ux_normed;
    
    dtxduxslip      = pslip.*(pystick.^2).*dtxduxstick./(uxynorm).^3.0;
    dtxduyslip      = -pslip.*(pxstick.*pystick).*dtyduystick./(uxynorm).^3.0;
    dtxdunslip      = mu.*kn(:,1).*kn(:,2).*un.^(kn(:,2)-1).*ux_normed;
    
    Ptx(stuck) 		= pxstick(stuck);
    Ptx(slipped) 	= pxslip(slipped);
    Ptx(sep) 		= 0.0;

    dtxdux(stuck)	= dtxduxstick(stuck);
    dtxdux(slipped)	= dtxduxslip(slipped);
    dtxdux(sep)		= 0.0;
    
    dtxduy(stuck)	= 0.0;
    dtxduy(slipped)	= dtxduyslip(slipped);
    dtxduy(sep)		= 0.0;

    dtxdun(stuck)	= dtxdunstick(stuck);
    dtxdun(slipped)	= dtxdunslip(slipped);
    dtxdun(sep)		= 0.0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Y-Tangential tractions & jacobians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pyslip = pslip.*uy_normed;
    
    dtyduxslip      = -pslip.*(pxstick.*pystick).*dtxduxstick./(uxynorm).^3.0;
    dtyduyslip      = pslip.*(pxstick.^2).*dtyduystick./(uxynorm).^3.0;
    dtydunslip      = mu.*kn(:,1).*kn(:,2).*un.^(kn(:,2)-1).*uy_normed;    

    Pty(stuck) 		= pystick(stuck);
    Pty(slipped) 	= pyslip(slipped);
    Pty(sep) 		= 0.0;
    
    dtydux(stuck)	= 0.0;
    dtydux(slipped)	= dtyduxslip(slipped);
    dtydux(sep)		= 0.0;
    
    dtyduy(stuck)	= dtyduystick(stuck);
    dtyduy(slipped)	= dtyduyslip(slipped);
    dtyduy(sep)		= 0.0;
    
    dtydun(stuck)	= dtydunstick(stuck);
    dtydun(slipped)	= dtydunslip(slipped);
    dtydun(sep)		= 0.0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     disp([num2str(length(sep)/Np*100) ' percentage separated'])
%     disp([num2str(length(slipped)/Np*100) ' percentage slipped'])
%     disp([num2str(length(stuck)/Np*100) ' percentage stuck'])
end