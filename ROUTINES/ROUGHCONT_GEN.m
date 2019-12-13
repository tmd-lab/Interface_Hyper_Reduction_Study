function [Ptx,Pty,Pn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun] ...
        = ROUGHCONT_GEN(uxyn,pn0,cpars,uxyntxynp) 
%JENKINS_GEN 3D Rough Contact model Pressure and Jacobians
% USAGE:
%	[Ptx,Pty,Pn,Jttx,Jtnx,Jtty,Jtny,Jnn] = ROUGHCONT_GEN(uxyn,pn0,pars,uxyntxynp);
% INPUTS:
%   uxyn	: (Npx3) 2 tangential and 1 normal displacement vectors 
%   pn0         : Static Normal pressure
%   pars        : Structure containing:
%       mu	: coefficient of friction (vector)
%       Tx	: x tangential stiffness properties (structure with a, b, c, d)
%       Ty	: y tangential stiffness properties (structure with a, b, c, d)
%       NF	: normal traction properties (structure with a, b, c, d, e)
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
    TX  = cpars.Tx;
    TY  = cpars.Ty;
    NF  = cpars.N;
    
    ux  = uxyn(:,1);
    uy  = uxyn(:,2);
    un  = uxyn(:,3) - cpars.Zn;
    uxp = uxyntxynp(:,1);
    uyp = uxyntxynp(:,2);
    unp = uxyntxynp(:,3) - cpars.Zn;
    
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
    pn = pn0 + exp(NF.a + NF.b.*log(un) + NF.c.*exp(-0.5*((log(un)+NF.d)./NF.e).^2));
    pn(sep) = 0;
    pslip   = mu.*pn;                               % Slip traction

%     unc = un;  unc(sep) = 0;  % un but zero where in separation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Stuck Predictors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% X %%%%%%%%%%%%%%%%%%
    ktx = exp(TX.a + TX.b.*tanh(TX.d.*(log(un)-TX.c)));
    dktxdn = TX.b.*TX.d.*(1.0-tanh(TX.d.*(log(un)-TX.c)).^2).*ktx./un;
    dktx2dn2 = -(2*TX.d./TX.b.*(log(ktx)-TX.a)+1).*dktxdn./un + (1./ktx).*dktxdn.^2;

    pxstick = ktx.*(ux-uxp) + dktxdn.*(un-unp).*ux;
%     pxstick(prevstick) = pxstick(prevstick) + txp(prevstick); % Memory preserved if it was stuck at turning point
    pxstick = pxstick + txp;  % CORRECTION: It's always stuck at the turning point because its TURNING!
    
    dtxduxstick = ktx + dktxdn.*(un-unp);
    dtxdunstick = dktxdn.*(2*ux-uxp) + dktx2dn2.*(un-unp).*ux;
    % Predict zeros for separated cases
    ktx(sep) = 0;       dktxdn(sep) = 0;        dktx2dn2(sep) = 0;
    pxstick(sep) = 0;   dtxduxstick(sep) = 0;   dtxdunstick(sep) = 0;
    
    %%%%%%%%%%%%%%%%%% Y %%%%%%%%%%%%%%%%%%
    kty = exp(TY.a + TY.b.*tanh(TY.d.*(log(un)-TY.c)));
    dktydn = TY.b.*TY.d.*(1.0-tanh(TY.d.*(log(un)-TY.c)).^2).*kty./un;
    dkty2dn2 = -(2*TY.d./TY.b.*(log(kty)-TY.a)+1).*dktydn./un + (1./kty).*dktydn.^2;

    pystick = kty.*(uy-uyp) + dktydn.*(un-unp).*uy;
%     pystick(prevstick) = pystick(prevstick) + typ(prevstick); % Memory preserved if it was stuck at turning point
    pystick = pystick + typ;  % CORRECTION: It's always stuck at the turning point because its TURNING!
    
    dtyduystick = kty + dktydn.*(un-unp);
    dtydunstick = dktydn.*(2*uy-uyp) + dkty2dn2.*(un-unp).*uy;
    % Predict zeros for separated cases
    kty(sep) = 0;       dktydn(sep) = 0;        dkty2dn2(sep) = 0;
    pystick(sep) = 0;   dtyduystick(sep) = 0;   dtydunstick(sep) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Setup Slip directions appropriately %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uxynorm = sqrt(pxstick.^2+pystick.^2);
    uxynorm(abs(uxynorm)<eps) = 1.0;
    ux_normed = pxstick./uxynorm;  ux_normed(abs(ux_normed)<eps) = 0.0;
    uy_normed = pystick./uxynorm;  uy_normed(abs(uy_normed)<eps) = 0.0;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine Stuck & Slipped Cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slipped = find(sqrt(pxstick.^2+pystick.^2)>abs(pslip));
    slipped = setdiff(slipped, sep);
    stuck = setdiff(setdiff(1:Np,slipped),sep);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Normal Traction & Jacobian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pn = pn - pn0;
    dtndun = (NF.b - NF.c./(NF.e.^2).*(log(un)+NF.d).*exp(-0.5*((log(un)+NF.d)./NF.e).^2)).*(pn./un);
    dtndun(sep) = 0.0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % X-Tangential tractions & jacobians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pxslip          = pslip.*ux_normed;

    dtxduxslip      = pslip.*(pystick.^2).*dtxduxstick./(uxynorm.^3.0);
    dtxduyslip      = -pslip.*(pxstick.*pystick).*dtyduystick./(uxynorm.^3.0);
    dtxdunslip      = mu.*dtndun.*ux_normed +...
        pslip.*(pystick.^2.*dtxdunstick - pxstick.*pystick.*dtydunstick)./(uxynorm.^3.0);

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
    pyslip          = pslip.*uy_normed;

    dtyduxslip      = -pslip.*(pxstick.*pystick).*dtxduxstick./(uxynorm.^3.0);
    dtyduyslip      = pslip.*(pxstick.^2).*dtyduystick./(uxynorm.^3.0);
    dtydunslip      = mu.*dtndun.*uy_normed +...
        pslip.*(pxstick.^2.*dtydunstick - pxstick.*pystick.*dtxdunstick)./(uxynorm.^3.0);

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