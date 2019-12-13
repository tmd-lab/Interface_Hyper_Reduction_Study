function [HYST,BACKBONE,UNLOAD,RELOAD] = GENSTEPPEDHYSTERESIS(nlresfun, contactfun, ...
                                             Xstart, Alphamax, Np, ...
                                             MESH, L, opt, QuadMats, varargin)
%STEPPEDHYSTERESIS Quasi-statically extracts the Hysteresis and
%Backbone
% USAGE:
%	[HYST,BACKBONE] = STEPPEDHYSTERESIS(nlresfun, nlforcefun,
%	Xstart, Alphamax, Np, MESH, L, No, opt, F2P);
% INPUTS:
%   nlresfun	:
%   nlforcefun	:
%   Xstart	:
%   Alphamax	:
%   Np		:
%   MESH	:
%   L		:
%   opt		:
%   F2P     :
% OUTPUTS:
%   HYST	:
%   BACKBONE	:

    Nqdofs = size(QuadMats.Q,1);
    Nq = Nqdofs/3;
    %% Backbone
    prev.uxyntxyn	= zeros(Nq,6);
    
    ALPHAS_BB	= linspace(0.0,Alphamax,Np);
    Na		= length(ALPHAS_BB);
    Ubb		= zeros(length(Xstart),Na);
    Pbb     = zeros(Nq*3,Na);
    Xguess 	= Xstart;
	if nargin == 10
        disp(' --------------------------------- ')
        disp(' Initializing Backbone Calculation ')
        disp(' --------------------------------- ')
		fprintf('Out of: %d\n', Na);
    end
    for i=1:Na
%         opt.FunctionTolerance = (ALPHAS_BB(i)+eps)^4;
        [Ubb(:,i),~,eflg]      = fsolve(@(X) nlresfun(X,ALPHAS_BB(i),prev), Xguess, ...
            opt);
        if eflg<=0
            fprintf('U');
            if i>1
                Xguess = Ubb(:,i-1);
            else
                Xguess = Xstart;
            end
            [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_BB(i),prev), Xguess, ...
                opt);
            if eflg<=0 && i>1
                Xguess = Xstart;
                [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_BB(i),prev), Xguess, ...
                    opt);
                if eflg<=0
                    error('What the hell.');
                end
            end
        end
        Ueph          = L*Ubb(:,i);
        [Ptx,Pty,Pn]  = contactfun(reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)',prev.uxyntxyn);
        Pbb(:,i)      = reshape([Ptx Pty Pn]',Nq*3,1);
        [~,dRdX,dRda] = nlresfun(Ubb(:,i),ALPHAS_BB(i),prev);
        dXda = -dRdX\dRda;
        if i<Na % First order predictor as guess
%             Xguess 	  = Ubb(:,i) + dXda*(ALPHAS_BB(i+1)-ALPHAS_BB(i));
            dX = dXda*(ALPHAS_BB(i+1)-ALPHAS_BB(i));
            [r, j] = nlresfun(Ubb(:, i), ALPHAS_BB(i+1), prev); e0 = r'*(-j\r);
            [r, j] = nlresfun(Ubb(:, i) + dX, ALPHAS_BB(i+1), prev); e1 = r'*(-j\r);
            if (-e0/(e1-e0))<1.01 && (-e0/(e1-e0))>-0.01
                Xguess = Ubb(:, i) - e0/(e0-e1)*dX;
            elseif abs(e0)<abs(e1)
                Xguess = Ubb(:,i);
            elseif abs(e1)<abs(e0)
                Xguess = Ubb(:,i) + dX;
            end
        end
        if nargin == 10
            % disp([num2str(i) '/' num2str(Na) ' DONE'])
			fprintf('%d; ', i);
        end
    end
    if nargin == 10
		fprintf('\n');
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
        disp(' Backbone Calculation Complete ')
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
    end
    
    %% Unloading
    prev.uxyntxyn(:,1:3) = reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)';
    prev.uxyntxyn(:,4:6) = [Ptx Pty Pn];
    
    ALPHAS_UL	= [ALPHAS_BB(end:-1:1) -ALPHAS_BB(2:end)];
    Naul	= length(ALPHAS_UL);
    Uul		= zeros(length(Xstart),Naul);
    Pul		= zeros(Nq*3,Naul);
    Xguess  = Ubb(:,end);
	if nargin == 10
        disp(' ---------------------------------- ')
        disp(' Initializing Unloading Calculation ')
        disp(' ---------------------------------- ')
		fprintf('Out of %d\n', Naul);
    end
    for i=1:Naul%
%         opt.FunctionTolerance = (ALPHAS_UL(i)+eps)^2;
        [Uul(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_UL(i),prev), Xguess, ...
            opt);
        if eflg<=0
            fprintf('U');
            if i>1
                Xguess = Uul(:,i-1);
            else
                Xguess = Ubb(:,end);
            end
            [Uul(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_UL(i),prev), Xguess, ...
                opt);
            if eflg<=0 && i>1
                Xguess = Ubb(:,end);
                [Uul(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_UL(i),prev), Xguess, ...
                    opt);
                if eflg<=0
                    error('Consider getting a life.');
                end
            end
        end
        Ueph 	 = L*Uul(:,i);
        [Ptx,Pty,Pn] = contactfun(reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)',prev.uxyntxyn);
        Pul(:,i)     = reshape([Ptx Pty Pn]',Nq*3,1);
        [~,dRdX,dRda] = nlresfun(Uul(:,i),ALPHAS_UL(i),prev);
        dXda = -dRdX\dRda;
        if i<Na % First order predictor as guess
%             Xguess 	 = Uul(:,i) + dXda*(ALPHAS_UL(i+1)-ALPHAS_UL(i));
            dX = dXda*(ALPHAS_UL(i+1)-ALPHAS_UL(i));
            [r, j] = nlresfun(Uul(:, i), ALPHAS_UL(i+1), prev); e0 = r'*(-j\r);
            [r, j] = nlresfun(Uul(:, i) + dX, ALPHAS_UL(i+1), prev); e1 = r'*(-j\r);
            if (-e0/(e1-e0))<1.01 && (-e0/(e1-e0))>-0.01
                Xguess = Ubb(:, i) - e0/(e0-e1)*dX;
            elseif abs(e0)<abs(e1)
                Xguess = Uul(:,i);
            elseif abs(e1)<abs(e0)
                Xguess = Uul(:,i) + dX;
            end
        end
        if nargin == 10
            % disp([num2str(i) '/' num2str(Naul) ' DONE'])
			fprintf('%d; ', i);
        end
    end
    if nargin == 10
		fprintf('\n');
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
        disp(' Unloading Calculation Complete ')
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
    end
    
    %% Reloading
    prev.uxyntxyn(:,1:3) = reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)';
    prev.uxyntxyn(:,4:6) = [Ptx Pty Pn];
    
    ALPHAS_RL	= ALPHAS_UL(end:-1:1);
    Narl	= length(ALPHAS_RL);
    Url		= zeros(length(Xstart),Narl);
    Prl		= zeros(Nq*3,Narl);
    Xguess  = Uul(:,end);
	if nargin == 10
        disp(' ---------------------------------- ')
        disp(' Initializing Reloading Calculation ')
        disp(' ---------------------------------- ')
		fprintf('Out of %d\n', Narl);
    end
    for i=1:Narl
%         opt.FunctionTolerance = (ALPHAS_RL(i)+eps)^2;
        [Url(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_RL(i),prev), Xguess, ...
                          opt);
        if eflg<=0
            fprintf('U');
            if i>1
                Xguess = Url(:,i-1);
            else
                Xguess = Uul(:,end);
            end
            [Url(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_RL(i),prev), Xguess, ...
                opt);
            if eflg<=0 && i>1
                Xguess = Uul(:,end);
                [Url(:,i),~,eflg] = fsolve(@(X) nlresfun(X,ALPHAS_RL(i),prev), Xguess, ...
                    opt);
                if eflg<=0
                    error('This is life.');
                end
            end
        end
        Ueph 	 = L*Url(:,i);
        [Ptx,Pty,Pn] = contactfun(reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)',prev.uxyntxyn);
        Prl(:,i)     = reshape([Ptx Pty Pn]',Nq*3,1);
        [~,dRdX,dRda] = nlresfun(Url(:,i),ALPHAS_RL(i),prev);
        dXda = -dRdX\dRda;
        if i<Na % First order predictor as guess
            dX = dXda*(ALPHAS_RL(i+1)-ALPHAS_RL(i));
            [r, j] = nlresfun(Url(:, i), ALPHAS_RL(i+1), prev); e0 = r'*(-j\r);
            [r, j] = nlresfun(Url(:, i) + dX, ALPHAS_RL(i+1), prev); e1 = r'*(-j\r);
            if (-e0/(e1-e0))<1.01 && (-e0/(e1-e0))>-0.01
                Xguess = Ubb(:, i) - e0/(e0-e1)*dX;
            elseif abs(e0)<abs(e1)
                Xguess = Url(:,i);
            elseif abs(e1)<abs(e0)
                Xguess = Url(:,i) + dX;
            end
        end
        if nargin == 10
            % disp([num2str(i) '/' num2str(Narl) ' DONE'])
			fprintf('%d; ', i);
        end
    end
    if nargin == 10
		fprintf('\n');
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
        disp(' Reloading Calculation Complete ')
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
    end
    %% Report Results
    BACKBONE.U	= Ubb;
    BACKBONE.P  = Pbb;
    BACKBONE.F	= QuadMats.T*Pbb;
    BACKBONE.A  = ALPHAS_BB;

    UNLOAD.U	= Uul;
    UNLOAD.P    = Pul;    
    UNLOAD.F	= QuadMats.T*Pul;
    UNLOAD.A    = ALPHAS_UL;
    
    RELOAD.U	= Url;
    RELOAD.P    = Prl;    
    RELOAD.F	= QuadMats.T*Prl;
    RELOAD.A    = ALPHAS_RL;
    
    HYST.U	= [(Uul(:,1)+Url(:,end))/2, Uul(:,2:end-1),...
                   (Uul(:,end)+Url(:,1))/2, Url(:,2:end-1)];
    HYST.P	= [(Pul(:,1)+Prl(:,end))/2, Pul(:,2:end-1),...
                   (Pul(:,end)+Prl(:,1))/2, Prl(:,2:end-1)];
    HYST.F  = QuadMats.T*HYST.P;               
    HYST.A  = [ALPHAS_UL ALPHAS_RL(2:end-1)];
    

end    