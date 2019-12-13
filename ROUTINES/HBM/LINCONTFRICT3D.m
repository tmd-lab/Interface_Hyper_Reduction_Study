function [FX,FY,FN,JXX,JXY,JXN,JYX,JYY,JYN,JNN] ...
        = LINCONTFRICT3D(UX,UY,UN,h,p0,pars,Nt)
%LINCONTFRICT3D Returns the forcing and jacobian for the 3D
%linearly compliant contact model with a penalty spring in the
%normal direction and coupled tangential friction model in the x &
%y directions 
% USAGE:
%	[FX,FY,FN,JXX,JXY,JXN,JYX,JYY,JYN,JNN] = ...
%		LINCONTFRICT3D(UX,UY,UN,h,f0,pars,Nt);
% INPUTS:
%   UX,UY,UN		: (NhcxNp) harmonic coefficients 
%   h			: (Nhcx1) vector of harmonic coefficients
%   f0			: (Npx1) static forces    
%   pars		: Structure of parameters    
%   Nt			: (int) Number of points for AFT
% OUTPUTS:
%   FX,FY,FN,JXX,JXY,JXN,JYX,JYY,JYN,JNN

    %%%%%%%%%%%%%%%%%%%
    % Interpretations %
    %%%%%%%%%%%%%%%%%%%
    Np	= size(UX,2);	% Number of contact points
    
    Nh	= length(h);	% Number of harmonics
    Nhc	= 0;		% Number of harmonic components
    zex	= 1;		% 1 if Zero harmonic exists; 0 if not
    if h(1)==0 		% Zero harmonic is present
        Nh 	= Nh-1;
        Nhc	= 2*Nh+1;
    else 		% Zero harmonic not present
        Nhc 	= 2*Nh;
        zex	= 0;
    end
    if size(UX,1)~=Nhc
        error('Input is Garbage.');
    end
    
    %%%%%%%%%%%%%%%%%%%
    % Initializations %
    %%%%%%%%%%%%%%%%%%%
    fx	= zeros(Nt,Np);		Fx	= zeros(Nhc,Np);
    fy	= zeros(Nt,Np);		Fy	= zeros(Nhc,Np);
    ftn = zeros(Nt,Np); 	Fn	= zeros(Nhc,Np); % forces
    
    jxx	= zeros(Nt,Nhc,Np);	JXX	= zeros(Nhc,Nhc,Np);
    jxy = zeros(Nt,Nhc,Np);	JXY	= zeros(Nhc,Nhc,Np);
    jxn	= zeros(Nt,Nhc,Np);	JXN	= zeros(Nhc,Nhc,Np); % x
                                                             % force jacobians
    jyx	= zeros(Nt,Nhc,Np);	JYX	= zeros(Nhc,Nhc,Np);
    jyy = zeros(Nt,Nhc,Np);	JYY	= zeros(Nhc,Nhc,Np);
    jyn	= zeros(Nt,Nhc,Np);	JYN	= zeros(Nhc,Nhc,Np); % y
                                                             % force jacobians
    jnn = zeros(Nt,Nhc,Np);	JNN	= zeros(Nhc,Nhc,Np); % normal force jacobians
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Frequency to Time Domain Transformation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h		= reshape(h,1,length(h));
    uxyn	= TIMESERIES_DERIV(Nt,h,[UX UY UN],0);
    ux		= uxyn(:,1:Np);
    uy		= uxyn(:,Np+(1:Np));
    un		= uxyn(:,2*Np+(1:Np));
    
    cst		= TIMESERIES_DERIV(Nt,h,eye(Nhc),0); % Fourier basis
                                                 % components
                                                 % in time
                                                 % domain
    cstmd	= cst.*ones(Nt,Nhc,Np);
    %%%%%%%%%%%%%%%%%%%%%%
    % Contact Parameters %
    %%%%%%%%%%%%%%%%%%%%%%
    mu	= pars.mu;	% (1x1) or (Npx1)
    kn	= pars.kn;	% (1x1) or (Npx1)
    ktx = pars.ktx;	% (1x1) or (Npx1)
    kty	= pars.kty;	% (1x1) or (Npx1)
                        % p0 also (1x1) or (Npx1)
    tol	= pars.tol;	% (1x1) always
    if length(mu)==1
        mu  = mu(1)*ones(Np,1);
        kn  = kn(1)*ones(Np,1);
        ktx = ktx(1)*ones(Np,1);
        kty = kty(1)*ones(Np,1);
    else
        if length(mu)~=Np || length(kn)~=Np || length(ktx)~=Np || ...
                length(kty)~=Np
            error('Garbage parameter inputs.');
        end
    end
    %%% REMOVE
    t=linspace(0,2*pi,Nt+1)';t(end)=[];
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Normal Contact Model %
    %%%%%%%%%%%%%%%%%%%%%%%%
    fn		= max(0,p0'+kn'.*un);
    fslip	= mu'.*fn; % move out
    sepcont	= (fn~=0); % 0 for separation; 1 for contact
    sepcontmd = permute(sepcont.*ones(Nt,Np,Nhc),[1 3 2]);
    
    jnn 	= (reshape(kn,1,1,Np).*cstmd).*sepcontmd;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tangential Contact Model %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stickslip	= false(Nt,Np);
    fxstick	= zeros(Nt,Np);	fystick		= zeros(Nt,Np);
    fxslip	= zeros(Nt,Np);	fyslip 		= zeros(Nt,Np);
    fnorm	= zeros(Nt,Np);
    fxnorm	= zeros(Nt,Np);	fynorm		= zeros(Nt,Np);
    fnormprev	= zeros(Nt,Np); fnormstick	= zeros(Nt,Np);
    
    jxxstick	= zeros(Nt,Nhc,Np);	jxystick	= zeros(Nt,Nhc,Np);
    jxxslip	= zeros(Nt,Nhc,Np);	jxyslip 	= zeros(Nt,Nhc,Np);
    jxnstick	= zeros(Nt,Nhc,Np);	jxnslip		= zeros(Nt,Nhc,Np);
    
    jyxstick	= zeros(Nt,Nhc,Np);	jyystick	= zeros(Nt,Nhc,Np);
    jyxslip	= zeros(Nt,Nhc,Np);	jyyslip 	= zeros(Nt,Nhc,Np);
    jynstick	= zeros(Nt,Nhc,Np);	jynslip		= zeros(Nt,Nhc,Np);
    m1		= @(l) mod((l-1)-1,Nt)+1;
    
    % Initialization for iterates
    fx      = ktx'.*(ux-UX(1,:));
    fy      = kty'.*(uy-UY(1,:));
    fnorm   = sqrt(fx.^2+fy.^2);
    fn0     = find(abs(fnorm)<eps);
    fxnorm  = fx./fnorm;    fxnorm(fn0)     = 0.0;
    fynorm  = fy./fnorm;    fynorm(fn0)     = 0.0;
    for i=1:Np
        jxx(:,1:end,i)     = ktx(i);
        jyy(:,1:end,i)     = kty(i);
    end
    it = 0;
    
    relerrp = 1.0;
    while 1 % Infinite loop for conducting time marching on
            % provided time series
        fnormprev	= fnorm;
        fxprev		= fx;
        fyprev		= fy;
        for n=1:Nt
            % FORCES
            fxstick(n,:)	= ktx'.*(ux(n,:)-ux(m1(n),:))+fx(m1(n),:);
            fystick(n,:)	= kty'.*(uy(n,:)-uy(m1(n),:))+fy(m1(n),:);
            fnorm(n,:)      = sqrt(fxstick(n,:).^2+fystick(n,:).^2); % Predicted fnorm
            fn0             = find(fnorm(n,:)==0);
            fxnorm(n,:)     = fxstick(n,:)./fnorm(n,:); fxnorm(n,fn0)   = 0.0;
            fynorm(n,:)     = fystick(n,:)./fnorm(n,:); fynorm(n,fn0)   = 0.0;
            
            fxslip(n,:)	= fslip(n,:).*fxnorm(n,:);
            fyslip(n,:)	= fslip(n,:).*fynorm(n,:);
            
            stickslip(n,:)	= (fnorm(n,:)>fslip(n,:)); % 0 for stick; 1 for
                                                       % slip
            fx(n,:)		= (fxstick(n,:).*(~stickslip(n,:)) + fxslip(n,:).*(stickslip(n,:))).*sepcont(n,:);
            fy(n,:)		= (fystick(n,:).*(~stickslip(n,:)) + fyslip(n,:).*(stickslip(n,:))).*sepcont(n,:);
            fnormstick(n,:) 	= fnorm(n,:);
            fnorm(n,:)      = sqrt(fx(n,:).^2+fy(n,:).^2); % Actual fnorml
            fn0             = find(fnorm(n,:)==0);
            fxnorm(n,:)     = fx(n,:)./fnorm(n,:); fxnorm(n,fn0)   = 0.0;
            fynorm(n,:)     = fy(n,:)./fnorm(n,:); fynorm(n,fn0)   = 0.0;
            
            % JACOBIANS
            jxxstick(n,:,:)     = reshape(ktx,[1 1 Np]).*(cstmd(n,:,:)-cstmd(m1(n),:,:))+jxx(m1(n),:,:);
            jxystick(n,:,:)     = jxy(m1(n),:,:);
            jxnstick(n,:,:)     = jxn(m1(n),:,:);
            
            jyxstick(n,:,:)     = jyx(m1(n),:,:);
            jyystick(n,:,:) 	= reshape(kty,[1 1 Np]).*(cstmd(n,:,:)-cstmd(m1(n),:,:))+jyy(m1(n),:,:);
            jynstick(n,:,:)     = jyn(m1(n),:,:);            
            
            jxxslip(n,:,:)      = reshape(fslip(n,:),[1 1 Np]).*( reshape(fystick(n,:).^2,[1 1 Np]).*jxxstick(n,:,:) - reshape(fxstick(n,:).*fystick(n,:),[1 1 Np]).*jyxstick(n,:,:) )./...
                reshape(fnormstick(n,:).^3,[1 1 Np]);
            jxyslip(n,:,:)      = reshape(fslip(n,:),[1 1 Np]).*( reshape(fystick(n,:).^2,[1 1 Np]).*jxystick(n,:,:) - reshape(fxstick(n,:).*fystick(n,:),[1 1 Np]).*jyystick(n,:,:) )./...
                reshape(fnormstick(n,:).^3,[1 1 Np]);
            jxnslip(n,:,:)      = reshape((mu.*kn)'.*fxnorm(n,:),[1 1 Np]).*cstmd(n,:,:) + ...
                reshape(fslip(n,:)./(fnormstick(n,:).^3),[1 1 Np]).*( reshape(fystick(n,:).^2,[1 1 Np]).*jxnstick(n,:,:) - reshape(fxstick(n,:).*fystick(n,:),[1 1 Np]).*jynstick(n,:,:) );
            fn0			= find(fnormstick(n,:)==0);
            jxxslip(n,:,fn0) = 0.0;
            jxyslip(n,:,fn0) = 0.0;
            jxnslip(n,:,fn0) = 0.0;
            
            jyxslip(n,:,:)      = reshape(fslip(n,:),[1 1 Np]).*( reshape(fxstick(n,:).^2,[1 1 Np]).*jyxstick(n,:,:) - reshape(fxstick(n,:).*fystick(n,:),[1 1 Np]).*jxxstick(n,:,:) )./...
                reshape(fnormstick(n,:).^3,[1 1 Np]);
            jyyslip(n,:,:)      = reshape(fslip(n,:),[1 1 Np]).*( reshape(fxstick(n,:).^2,[1 1 Np]).*jyystick(n,:,:) - reshape(fxstick(n,:).*fystick(n,:),[1 1 Np]).*jxystick(n,:,:) )./...
                reshape(fnormstick(n,:).^3,[1 1 Np]);
            jynslip(n,:,:)      = reshape((mu.*kn)'.*fynorm(n,:),[1 1 Np]).*cstmd(n,:,:) + ...
                reshape(fslip(n,:)./(fnormstick(n,:).^3),[1 1 Np]).*( reshape(fxstick(n,:).^2,[1 1 Np]).*jynstick(n,:,:) - reshape(fxstick(n,:).*fystick(n,:),[1 1 Np]).*jxnstick(n,:,:) );
            jyxslip(n,:,fn0) = 0.0;
            jyyslip(n,:,fn0) = 0.0;
            jynslip(n,:,fn0) = 0.0;
            
            jxx(n,:,:)		= (jxxstick(n,:,:).*reshape(~stickslip(n,:),[1 1 Np]) + jxxslip(n,:,:).*reshape(stickslip(n,:),[1 1 Np])).*sepcontmd(n,:,:);
            jxy(n,:,:)		= (jxystick(n,:,:).*reshape(~stickslip(n,:),[1 1 Np]) + jxyslip(n,:,:).*reshape(stickslip(n,:),[1 1 Np])).*sepcontmd(n,:,:);
            jxn(n,:,:)		= (jxnstick(n,:,:).*reshape(~stickslip(n,:),[1 1 Np]) + jxnslip(n,:,:).*reshape(stickslip(n,:),[1 1 Np])).*sepcontmd(n,:,:);            
            
            jyx(n,:,:)		= (jyxstick(n,:,:).*reshape(~stickslip(n,:),[1 1 Np]) + jyxslip(n,:,:).*reshape(stickslip(n,:),[1 1 Np])).*sepcontmd(n,:,:);
            jyy(n,:,:)		= (jyystick(n,:,:).*reshape(~stickslip(n,:),[1 1 Np]) + jyyslip(n,:,:).*reshape(stickslip(n,:),[1 1 Np])).*sepcontmd(n,:,:);
            jyn(n,:,:)		= (jynstick(n,:,:).*reshape(~stickslip(n,:),[1 1 Np]) + jynslip(n,:,:).*reshape(stickslip(n,:),[1 1 Np])).*sepcontmd(n,:,:);
        end
        it = it+1;
%         upd8i = find(rms(fnorm-fnormprev)>tol | (rms(fnorm)>eps &
%         rms(fnorm-fnormprev)./rms(fnorm)>tol) ); % Selective update

        abserr = rms(fnorm(end,:)-fnormprev(end,:));
        relerr = abs(abserr/max(fnorm(end,:)));
        
%         abserr = max(abs(fnorm-fnormprev));
%         maf    = max(max(abs(fnorm)), eps);
%         relerr = abserr./maf;
        
%         fprintf('%d %e\n',it,max(relerr));
        if relerr<tol
            break;
        end
        
        if it==2
            break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time to Frequency Domain Transformation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FX	= GETFOURIERCOEFF(h,fx);
    FY	= GETFOURIERCOEFF(h,fy);
    FN	= GETFOURIERCOEFF(h,fn-p0');
    
    JXX	= reshape(GETFOURIERCOEFF(h,reshape(jxx,[Nt Nhc*Np])),[Nhc Nhc Np]);
    JXY	= reshape(GETFOURIERCOEFF(h,reshape(jxy,[Nt Nhc*Np])),[Nhc Nhc Np]);
    JXN	= reshape(GETFOURIERCOEFF(h,reshape(jxn,[Nt Nhc*Np])),[Nhc Nhc Np]);
    
    JYX	= reshape(GETFOURIERCOEFF(h,reshape(jyx,[Nt Nhc*Np])),[Nhc Nhc Np]);
    JYY	= reshape(GETFOURIERCOEFF(h,reshape(jyy,[Nt Nhc*Np])),[Nhc Nhc Np]);
    JYN	= reshape(GETFOURIERCOEFF(h,reshape(jyn,[Nt Nhc*Np])),[Nhc Nhc Np]);
    
    JNN	= reshape(GETFOURIERCOEFF(h,reshape(jnn,[Nt Nhc*Np])),[Nhc Nhc Np]);
    
    % Accounting for zero harmonic uncertainty for fully stuck case
    if h(1)==0
        for i=1:Np
            if isempty(find(stickslip(:,i)==1, 1))  % No sliding
                %#1. SET ALL TO 0
%                 JXX(1, :, i) = 0;
%                 JXY(1, :, i) = 0;
%                 JXN(1, :, i) = 0;
% 
%                 JYX(1, :, i) = 0;
%                 JYY(1, :, i) = 0;
%                 JYN(1, :, i) = 0;
%                 
%                 FX(1, i) = 0;
%                 JXX(1, 1, i) = ktx;
%                 
%                 FY(1, i) = 0;
%                 JYY(1, 1, i) = kty;

%                 %#2. MAKE STUCK CASE TO LINEAR SPRING
%                 JXX(:, :, i) = eye(Nhc)*ktx;
%                 JXY(:, :, i) = zeros(Nhc);
%                 JXN(:, :, i) = zeros(Nhc);
%                 FX = JXX*UX;
%                 
%                 JYX(:, :, i) = zeros(Nhc);
%                 JYY(:, :, i) = eye(Nhc)*kty;
%                 JXN(:, :, i) = zeros(Nhc);
%                 FY = JYY*UY;

%                 %#3. ADD LINEAR SPRING ON ZERO HARMONIC
%                 FX(1, i) = FX(1, i) + ktx*UX(1, i);
%                 JXX(1, 1, i) = ktx;
%                 
%                 FY(1, i) = FY(1, i) + kty*UY(1, i);
%                 JYY(1, 1, i) = kty;
            end
        end
    end
end