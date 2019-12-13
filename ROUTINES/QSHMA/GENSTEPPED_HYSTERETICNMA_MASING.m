function [BACKBONE] = GENSTEPPED_HYSTERETICNMA_MASING(nlresfun, contactfun, ...
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
    DS_BB   = zeros(size(QS_BB));
    Na		= length(QS_BB);
    Ubb		= zeros(length(Xstat)+1,Na);
    Pbb     = zeros(Nq*3,Na);
    Xstart = Xstat;
	Xguess 	= [Xstat+msp(1:end-1)*QS_BB(2); msp(end)]./Cmat(1:end-1);
    
    % LGL Quadrature Points
    [xi, wi] = LGLWT(Copt.Nqp, 0, 1);
    qr = zeros(size(xi));
    ar = zeros(size(xi));
    
    for i=2:Na
        DS_BB(i) = 0;
        for k=2:(Copt.Nqp+1)
            qr(k) = QS_BB(i)*xi(k);
            % LINEAR PREDICTOR
            Xguessp = [Xstat+msp(1:end-1)*qr(k); msp(end)]./Cmat(1:end-1);
            [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, qr(k)/Cmat(end), Xstart, prev), Xguessp, opt);
            Ubb(:,i) = Cmat(1:end-1).*Ubb(:,i);
            if eflg<=0  % LINE SEARCH PREDICTOR
                [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, qr(k)/Cmat(end), Xstart, prev), Xguess, opt);
                Ubb(:,i) = Cmat(1:end-1).*Ubb(:,i);
            end
            if eflg<=0  % FAILURE TO CONVERGE
                fprintf('U');
                if k>1  
                    Xguess = Ubb(:,k-1);
                else
                    Xguess = [Xstat+msp(1:end-1)*qr(k); msp(end)]./Cmat(1:end-1);
                end
                [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, qr(k)/Cmat(end), Xstart, prev), Xguess, opt);
                Ubb(:,i) = Cmat(1:end-1).*Ubb(:,i);

                if eflg<=0 && i>1
                    Xguess = [Xstat+msp(1:end-1)*qr(k); msp(end)]./Cmat(1:end-1);
                    [Ubb(:,i),~,eflg] = fsolve(@(X) nlresfun(X, qr(k)/Cmat(end), Xstart, prev), Xguess, opt);
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
            [~,dRdX,dRda] = nlresfun(Ubb(:,i)./Cmat(1:end-1), qr(k)/Cmat(end), Xstart, prev);
            dXda = -dRdX\dRda;
            if k<Copt.Nqp % First order predictor as guess
    %             Xguess 	  = Ubb(:,i) + dXda*(ALPHAS_BB(i+1)-ALPHAS_BB(i));
                dX = dXda*(qr(k+1)-qr(k))/Cmat(end);
                if Copt.lsrch==1
                    [r, j] = nlresfun(Ubb(:, i)./Cmat(1:end-1), qr(k+1)/Cmat(end), Xstart, prev); e0 = r'*(-j\r);
                    [r, j] = nlresfun(Ubb(:, i)./Cmat(1:end-1) + dX, qr(k+1)/Cmat(end), Xstart, prev); e1 = r'*(-j\r);
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
            
            DS_BB(i) = DS_BB(i) + QS_BB(i)*wi(k)*qr(k)*Ubb(end,i);
        end
        
        if nargin == 13
            % disp([num2str(i) '/' num2str(Na) ' DONE'])
			fprintf('%d; ', i);
        end
    end
 
    BACKBONE = 0;
end