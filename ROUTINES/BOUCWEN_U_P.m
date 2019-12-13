function [Ptx,Pty,Pn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun] ...
        = BOUCWEN_U_P(uxyn,pn0,pars,uxyntxynp) 
%BOUCWEN_U_P 3D Planar uncoupled Bouc-Wen model Pressure and Jacobians
% USAGE:
%	[Ptx,Pty,Pn,Jttx,Jtnx,Jtty,Jtny,Jnn] = BOUCWEN_U_P(uxyn,pn0,pars,uxyntxynp);
% INPUTS:
%   uxyn	: (Npx3) 2 tangential and 1 normal displacement vectors 
%   pn0         : Static Normal pressure
%   pars        : Structure containing:
%       mu	: coefficient of friction
%       sig	: parameter sigma
%       rho : parameter rho
%       ne  : parameter n (exponent power)
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
    sig = pars.sig;
    rho = pars.rho;
    ne  = pars.ne;

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
    pslip   = mu.*pn;                             % Slip traction
    pslipp  = mu.*(tnp+pn0);    % previous slip traction
    
    sep = find(pn==0); 		% Separated cases
    cnt = setdiff(1:Np, sep)'; 	% Contacted cases    
        
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
   
%     % Fully-Consistent Euler Step
% 	% X direction traction & Jacobians
%     txbypsl_p = txp./pslipp;  txbypsl_p(txp==0) = 0;
% 	tybypsl_p = typ./pslipp;  tybypsl_p(typ==0) = 0;
%     Ptx             = pslipp.*rho.*(1-(sig.*sign((ux-uxp).*(txp+eps))+(1-sig)).*(txbypsl_p).^ne).*(ux-uxp) + txp;
% 	Ptx(sep) 		= 0.0;
% 
% 	dtxdux          = pslipp.*rho.*(1-(sig.*sign((ux-uxp).*(txp+eps))+(1-sig)).*(txbypsl_p).^ne);
%     dtxdux(sep)		= 0.0;
% 
% 	dtxduy      	= 0.0;
% 
%     dtxdun          = 0.0;
% 
%     % Y direction traction & Jacobian
% 	Pty      		= pslipp.*rho.*(1-(sig.*sign((uy-uyp).*(typ+eps))+(1-sig)).*(tybypsl_p).^ne).*(uy-uyp) + typ;
%     Pty(sep) 		= 0.0;
% 
%     dtydux      	= 0.0;
% 
%     dtyduy          = pslipp.*rho.*(1-(sig.*sign((uy-uyp).*(typ+eps))+(1-sig)).*(tybypsl_p).^ne);
% 	dtyduy(sep)		= 0.0;
% 
%     dtydun      	= 0.0;
    
% 	% Half-Consistent Euler Step
% 	% X direction traction & Jacobians
%     txbypsl_p = txp./pslipp;  txbypsl_p(txp==0) = 0;
% 	tybypsl_p = typ./pslipp;  tybypsl_p(typ==0) = 0;
%     Ptx             = pslip.*rho.*(1-(sig.*sign((ux-uxp).*(txp+eps))+(1-sig)).*(txbypsl_p).^ne).*(ux-uxp) + txp;
% 	Ptx(sep) 		= 0.0;
% 
% 	dtxdux          = pslip.*rho.*(1-(sig.*sign((ux-uxp).*(txp+eps))+(1-sig)).*(txbypsl_p).^ne);
%     dtxdux(sep)		= 0.0;
% 
% 	dtxduy      	= 0.0;
% 
%     dtxdun          = (Ptx-txp)./un;
%     dtxdun(sep)     = 0.0;
% 
%     % Y direction traction & Jacobian
% 	Pty      		= pslip.*rho.*(1-(sig.*sign((uy-uyp).*(typ+eps))+(1-sig)).*(tybypsl_p).^ne).*(uy-uyp) + typ;
%     Pty(sep) 		= 0.0;
% 
%     dtydux      	= 0.0;
% 
%     dtyduy          = pslip.*rho.*(1-(sig.*sign((uy-uyp).*(typ+eps))+(1-sig)).*(tybypsl_p).^ne);
% 	dtyduy(sep)		= 0.0;
% 
%     dtydun      	= (Pty-typ)./un;
% 	dtydun(sep)		= 0.0;

	% Inconsistent formula 
    % X direction traction & Jacobians
	Ptx             = pslip.*rho.*(1-(sig.*sign((ux-uxp).*(txp+eps))+(1-sig)).*(txp./pslip).^ne).*(ux-uxp) + txp;
    Ptx(sep) 		= 0.0;

	dtxdux          = pslip.*rho.*(1-(sig.*sign((ux-uxp).*(txp+eps))+(1-sig)).*(txp./pslip).^ne);
    dtxdux(sep)		= 0.0;

	dtxduy      	= 0.0;

    dtxdun          = rho.*(mu.*kn+(sig.*sign((ux-uxp).*(txp+eps))+(1-sig)).*(ne-mu.*kn)).*(ux-uxp);
	dtxdun(sep)     = 0.0;

    % Y direction traction & Jacobian
	Pty      		= pslip.*rho.*(1-(sig.*sign((uy-uyp).*(typ+eps))+(1-sig)).*(typ./pslip).^ne).*(uy-uyp) + typ;
	Pty(sep) 		= 0.0;

    dtydux      	= 0.0;

	dtyduy          = pslip.*rho.*(1-(sig.*sign((uy-uyp).*(typ+eps))+(1-sig)).*(typ./pslip).^ne);
    dtyduy(sep)		= 0.0;

	dtydun      	= rho.*(mu.*kn+(sig.*sign((uy-uyp).*(typ+eps))+(1-sig)).*(ne-mu.*kn)).*(uy-uyp);
    dtydun(sep)		= 0.0;

%     disp([num2str(length(sep)/Np*100) ' percentage separated'])
%     disp([num2str(length(slipped)/Np*100) ' percentage slipped'])
%     disp([num2str(length(stuck)/Np*100) ' percentage stuck'])
end
