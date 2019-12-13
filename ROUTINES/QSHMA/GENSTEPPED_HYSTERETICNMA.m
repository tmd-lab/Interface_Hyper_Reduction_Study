function [HYST,BACKBONE,UNLOAD,RELOAD] = GENSTEPPED_HYSTERETICNMA(nlresfun, contactfun, ...
                                             Xstat, msp, Cmat, Qsmax, Np, ...
                                             MESH, L, opt, QuadMats, Copt, varargin)
%GENSTEPPEDHYSTERESIS_ADMS Quasi-statically extracts the Hysteresis and
%Backbone
% USAGE:
%	[HYST,BACKBONE] = GENSTEPPEDHYSTERESIS_ADMS(nlresfun, nlforcefun,
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
    
    QS_BB	= linspace(0.0,Qsmax,Np);
    Na		= length(QS_BB);
    Ubb		= zeros(length(Xstat)+1,Na);
    Pbb     = zeros(Nq*3,Na);
	if nargin == 13
        disp(' --------------------------------- ')
        disp(' Initializing Backbone Calculation ')
        disp(' --------------------------------- ')
		fprintf('Out of: %d\n', Na);
    end
    Xstart = Xstat;
	Xguess 	= [Xstat+msp(1:end-1)*QS_BB(2); msp(end)]./Cmat(1:end-1);
    for i=1:Na
        % FIRST TRY TO USE LINEAR PREDICTOR
        Xguessp 	= [Xstat+msp(1:end-1)*QS_BB(i); msp(end)]./Cmat(1:end-1);
        [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_BB(i)/Cmat(end), Xstart, prev), Xguessp, opt);
        Ubb(:,i) = Cmat(1:end-1).*Ubb(:,i);
        if eflg<=0  % LINE SEARCH PREDICTOR
            [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_BB(i)/Cmat(end), Xstart, prev), Xguess, opt);
            Ubb(:,i) = Cmat(1:end-1).*Ubb(:,i);
        end
        
        if eflg<=0  % FAILURE TO CONVERGE
            fprintf('U');
            if i>1  
                Xguess = Ubb(:,i-1);
            else
                Xguess = [Xstat+msp(1:end-1)*QS_BB(i); msp(end)]./Cmat(1:end-1);
            end
            [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_BB(i)/Cmat(end), Xstart, prev), Xguess, opt);
            Ubb(:,i) = Cmat(1:end-1).*Ubb(:,i);
            
            if eflg<=0 && i>1
                Xguess = [Xstat+msp(1:end-1)*QS_BB(i); msp(end)]./Cmat(1:end-1);
                [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_BB(i)/Cmat(end), Xstart, prev), Xguess, opt);
                Ubb(:,i) = Cmat(1:end-1)*Ubb(:,i);
                
                if eflg<=0
                    error('What the hell.');
                end
            end
        end
        
        % CONVERGENCE
        Ueph          = L*Ubb(1:end-1,i);
        [Ptx,Pty,Pn]  = contactfun(reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)',prev.uxyntxyn);
        Pbb(:,i)      = reshape([Ptx Pty Pn]',Nq*3,1);
        [~,dRdX,dRda] = nlresfun(Ubb(:,i)./Cmat(1:end-1), QS_BB(i)/Cmat(end), Xstart, prev);
        dXda = -dRdX\dRda;
        if i<Na % First order predictor as guess
%             Xguess 	  = Ubb(:,i) + dXda*(ALPHAS_BB(i+1)-ALPHAS_BB(i));
            dX = dXda*(QS_BB(i+1)-QS_BB(i))/Cmat(end);
            if Copt.lsrch==1
                [r, j] = nlresfun(Ubb(:, i)./Cmat(1:end-1), QS_BB(i+1)/Cmat(end), Xstart, prev); e0 = r'*(-j\r);
                [r, j] = nlresfun(Ubb(:, i)./Cmat(1:end-1) + dX, QS_BB(i+1)/Cmat(end), Xstart, prev); e1 = r'*(-j\r);
                if (-e0/(e1-e0))<1.01 && (-e0/(e1-e0))>-0.01
                    Xguess = Ubb(:, i) - e0/(e0-e1)*dX;
                elseif abs(e0)<abs(e1)
                    Xguess = Ubb(:,i);
                elseif abs(e1)<=abs(e0)
                    Xguess = Ubb(:,i) + dX;
                end
            else
                Xguess = Ubb(:,i) + dX;
            end
        end
        
        if nargin == 13
            % disp([num2str(i) '/' num2str(Na) ' DONE'])
			fprintf('%d; ', i);
        end
    end
    if nargin == 13
		fprintf('\n');
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
        disp(' Backbone Calculation Complete ')
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
    end
    
    %% Unloading
    prev.uxyntxyn(:,1:3) = reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)';
    prev.uxyntxyn(:,4:6) = [Ptx Pty Pn];
    
    QS_UL	= [QS_BB(end:-1:2) -QS_BB(2:end)];
    Naul	= length(QS_UL);
    Uul		= zeros(length(Xstat)+1, Naul);
    Pul		= zeros(Nq*3,Naul);
    Xguess  = Ubb(:,end);
	if nargin == 13
        disp(' ---------------------------------- ')
        disp(' Initializing Unloading Calculation ')
        disp(' ---------------------------------- ')
		fprintf('Out of %d\n', Naul);
    end
    Xstart = Xstat;
%     Xstart = Ubb(1:end-1,end);
    for i=1:Naul%
        % FIRST TRY TO USE LINEAR PREDICTOR
        Xguessp 	= [Xstat+msp(1:end-1)*QS_UL(i); msp(end)]./Cmat(1:end-1);
        [Uul(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_UL(i)/Cmat(end), Xstart, prev), Xguessp, opt);
        Uul(:,i) = Cmat(1:end-1).*Uul(:,i);
        if eflg<=0 % LINE SEARCH PREDICTOR
            [Uul(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_UL(i)/Cmat(end), Xstart, prev), Xguess, opt);
            Uul(:,i) = Cmat(1:end-1).*Uul(:,i);
        end
        
        if eflg<=0  % FAILURE TO CONVERGE
            fprintf('U');
            if i>1
                Xguess = Uul(:,i-1);
            else
                Xguess = Ubb(:,end);
            end
            [Uul(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_UL(i)/Cmat(end), Xstart, prev), Xguess, opt);
            Uul(:,i) = Cmat(1:end-1).*Uul(:,i);
        
            if eflg<=0 && i>1
                Xguess = Ubb(:,end);
                [Uul(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_UL(i)/Cmat(end), Xstart, prev), Xguess, opt);
                Uul(:,i) = Cmat(1:end-1).*Uul(:,i);
                
                if eflg<=0
                    error('Consider getting a life.');
                end
            end
        end
        
        % CONVERGENCE
        Ueph 	 = L*Uul(1:end-1,i);
        [Ptx,Pty,Pn] = contactfun(reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)',prev.uxyntxyn);
        Pul(:,i)     = reshape([Ptx Pty Pn]',Nq*3,1);
        [~,dRdX,dRda] = nlresfun(Uul(:,i)./Cmat(1:end-1),QS_UL(i)/Cmat(end),Xstart, prev);
        dXda = -dRdX\dRda;
        if i<Na % First order predictor as guess
%             Xguess 	 = Uul(:,i) + dXda*(ALPHAS_UL(i+1)-ALPHAS_UL(i));
            dX = dXda*(QS_UL(i+1)-QS_UL(i))/Cmat(end);
            if Copt.lsrch==1
                [r, j] = nlresfun(Uul(:, i)./Cmat(1:end-1), QS_UL(i+1)/Cmat(end), Xstart, prev); e0 = r'*(-j\r);
                [r, j] = nlresfun(Uul(:, i)./Cmat(1:end-1) + dX, QS_UL(i+1)/Cmat(end), Xstart, prev); e1 = r'*(-j\r);
                if (-e0/(e1-e0))<1.01 && (-e0/(e1-e0))>-0.01
                    Xguess = Ubb(:, i) - e0/(e0-e1)*dX;
                elseif abs(e0)<abs(e1)
                    Xguess = Uul(:,i);
                elseif abs(e1)<=abs(e0)
                    Xguess = Uul(:,i) + dX;
                end
            else
                Xguess = Uul(:,i) + dX;
            end
        end
        
        if nargin == 13
            % disp([num2str(i) '/' num2str(Naul) ' DONE'])
			fprintf('%d; ', i);
        end
    end
    if nargin == 13
		fprintf('\n');
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
        disp(' Unloading Calculation Complete ')
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
    end
    
    %% Reloading
    prev.uxyntxyn(:,1:3) = reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)';
    prev.uxyntxyn(:,4:6) = [Ptx Pty Pn];
    
    QS_RL	= QS_UL(end:-1:1);
    Narl	= length(QS_RL);
    Url		= zeros(length(Xstat)+1, Narl);
    Prl		= zeros(Nq*3,Narl);
    Xguess  = Uul(:,end);
	if nargin == 13
        disp(' ---------------------------------- ')
        disp(' Initializing Reloading Calculation ')
        disp(' ---------------------------------- ')
		fprintf('Out of %d\n', Narl);
    end
    Xstart = Xstat;
%     Xstart = Uul(1:end-1,end);
    for i=1:Narl
        % FIRST TRY TO USE LINEAR PREDICTOR
        Xguessp 	= [Xstat+msp(1:end-1)*QS_RL(i); msp(end)]./Cmat(1:end-1);
        [Url(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_RL(i)/Cmat(end), Xstart, prev), Xguessp, opt);
        Url(:,i) = Cmat(1:end-1).*Url(:,i);
        if eflg<=0 % LINE SEARCH PREDICTOR
            [Url(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_RL(i)/Cmat(end), Xstart, prev), Xguess, opt);
            Url(:,i) = Cmat(1:end-1).*Url(:,i);
        end
        
        if eflg<=0   % FAILURE TO CONVERGE
            fprintf('U');
            if i>1
                Xguess = Url(:,i-1);
            else
                Xguess = Uul(:,end);
            end
            [Url(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_RL(i)/Cmat(end), Xstart, prev), Xguess, opt);
            Url(:,i) = Cmat(1:end-1).*Url(:,i);

            if eflg<=0 && i>1
                Xguess = Uul(:,end);
                [Url(:,i),~,eflg] = fsolve(@(X) nlresfun(X, QS_RL(i)/Cmat(end), Xstart, prev), Xguess, opt);
                Url(:,i) = Cmat(1:end-1).*Url(:,i);
                if eflg<=0
                    error('This is life.');
                end
            end
        end
        
        % CONVERGENCE
        Ueph 	 = L*Url(1:end-1,i);
        [Ptx,Pty,Pn] = contactfun(reshape(QuadMats.Q*Ueph(1:(MESH.dpn*MESH.Nn)),3,Nq)',prev.uxyntxyn);
        Prl(:,i)     = reshape([Ptx Pty Pn]',Nq*3,1);
        [~,dRdX,dRda] = nlresfun(Url(:,i)./Cmat(1:end-1),QS_RL(i)/Cmat(end),Xstart, prev);
        dXda = -dRdX\dRda;
        if i<Na % First order predictor as guess
            dX = dXda*(QS_RL(i+1)-QS_RL(i));
            if Copt.lsrch==1
                [r, j] = nlresfun(Url(:, i)./Cmat(1:end-1), QS_RL(i+1)/Cmat(end), Xstart, prev); e0 = r'*(-j\r);
                [r, j] = nlresfun(Url(:, i)./Cmat(1:end-1) + dX, QS_RL(i+1)/Cmat(end), Xstart, prev); e1 = r'*(-j\r);
                if (-e0/(e1-e0))<1.01 && (-e0/(e1-e0))>-0.01
                    Xguess = Ubb(:, i) - e0/(e0-e1)*dX;
                elseif abs(e0)<abs(e1)
                    Xguess = Url(:,i);
                elseif abs(e1)<=abs(e0)
                    Xguess = Url(:,i) + dX;
                end
            else
                Xguess = Url(:,i) + dX;
            end
        end
        
        if nargin == 13
            % disp([num2str(i) '/' num2str(Narl) ' DONE'])
			fprintf('%d; ', i);
        end
    end
    if nargin == 13
		fprintf('\n');
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
        disp(' Reloading Calculation Complete ')
        disp(' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ')
    end
    %% Report Results
    BACKBONE.U	= Ubb;
    BACKBONE.P  = Pbb;
    BACKBONE.F	= QuadMats.T*Pbb;
    BACKBONE.Q  = QS_BB;
    BACKBONE.W  = sqrt(Ubb(end,:));
    BACKBONE.A  = BACKBONE.Q.*(BACKBONE.W.^2);
    BACKBONE.MS = (Ubb(1:end-1,:)-Xstat)./QS_BB;
    BACKBONE.MS(:,1) = msp(1:end-1);

    UNLOAD.U	= Uul;
    UNLOAD.P    = Pul;    
    UNLOAD.F	= QuadMats.T*Pul;
    UNLOAD.Q    = QS_UL;
    UNLOAD.W    = sqrt(Uul(end,:));
    UNLOAD.A    = UNLOAD.Q.*(UNLOAD.W.^2);
    UNLOAD.MS   = (Uul(1:end-1,:)-Xstat)./QS_UL;
    
    RELOAD.U	= Url;
    RELOAD.P    = Prl;    
    RELOAD.F	= QuadMats.T*Prl;
    RELOAD.Q    = QS_RL;
    RELOAD.W    = sqrt(Url(end,:));
    RELOAD.A    = RELOAD.Q.*(RELOAD.W.^2);
    RELOAD.MS   = (Url(1:end-1,:)-Xstat)./QS_RL;
    
    HYST.U	= [(Uul(:,1)+Url(:,end))/2, Uul(:,2:end-1),...
                   (Uul(:,end)+Url(:,1))/2, Url(:,2:end-1)];
    HYST.P	= [(Pul(:,1)+Prl(:,end))/2, Pul(:,2:end-1),...
                   (Pul(:,end)+Prl(:,1))/2, Prl(:,2:end-1)];
    HYST.F  = QuadMats.T*HYST.P;               
    HYST.Q  = [QS_UL QS_RL(2:end-1)];
    HYST.W  = sqrt(HYST.U(end,:));
    HYST.A  = HYST.Q.*(HYST.W.^2);
    HYST.MS = [(UNLOAD.MS(:,1)+RELOAD.MS(:,end))/2, UNLOAD.MS(:,2:end-1),...
                   (UNLOAD.MS(:,end)+RELOAD.MS(:,1))/2, RELOAD.MS(:,2:end-1)];
    

end