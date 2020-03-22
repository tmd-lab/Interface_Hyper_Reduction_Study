function [ErrsWZD, dErrsdp, BB] = WJMODEL_BBFUN_LOOP(pars, mdi, expdat, K, M, X0, Fv, L, QuadMats, CFUN, Npatches, Nqp, opt)
%WJMODEL_BBFUN_LOOP returns the backbone by doing a full hysteresis loop to
%obtain damping. 
%
%Uses quadrature to integrate along the hysteresis loop

	% Prestress - From WJ with modifications
    
    prevS.uxynFxyn = zeros(Npatches,6);
    prevS.duxdp = zeros(Npatches,numel(pars));
    prevS.duydp = zeros(Npatches,numel(pars));
    prevS.dundp = zeros(Npatches,numel(pars));
    prevS.dFxyn = zeros(Npatches,3, 3);
    prevS.dFxyndpar = zeros(Npatches,numel(pars), 3);
    prevS.uxykp0 = zeros(size(prevS.uxynFxyn, 1), 2); %Tracks reference point of kp spring
    
    
    %CFUN for static analysis
    CFUNS = @(uxyn, pars) CFUN(uxyn, pars, prevS);
    
    [Xstat, eflg] = fsolve(@(X) WJRES_STR([X; 0.0], pars, K, 0*Fv, Fv, L, QuadMats.T, QuadMats.Q, CFUNS, Npatches), X0, opt);
    if eflg<0
        disp(pars)
        error('Non-convergent!');
    end
    [~,dRstat,~,dRdpstat] = WJRES_STR([Xstat; 0.0], pars, K, 0*Fv, Fv, L, QuadMats.T, QuadMats.Q, CFUNS, Npatches);
    dXdpstat = -dRstat\dRdpstat;
    
    uxyn = [QuadMats.Q(1:3:end,:)*Xstat, QuadMats.Q(2:3:end,:)*Xstat, QuadMats.Q(3:3:end,:)*Xstat];
    
    [Fxyn, prevS.dFxyn, prevS.dFxyndpar] = CFUNS(uxyn, pars);
	prevS.uxynFxyn = [uxyn, Fxyn];

% 	if size(dFxynstat,3)==1
%        tmp = dFxynstat;
%        dFxynstat = zeros(size(dFxynstat,1),size(dFxynstat,2),3);
%        dFxynstat(:,1,1) = tmp(:,1);
%        dFxynstat(:,2,2) = tmp(:,2);
%        dFxynstat(:,3,3) = tmp(:,3); % Np x uxyn x fxyn
%     end
% 	dFxynstat = num2cell(permute(dFxynstat,[3,2,1]), [2,1]);
%     prevS.dFxyndXdpar = blkdiag(dFxynstat{:})*QuadMats.Q*dXdpstat;
    
    prevS.duxdp = QuadMats.Q(1:3:end,:)*dXdpstat;
    prevS.duydp = QuadMats.Q(2:3:end,:)*dXdpstat;
    prevS.dundp = QuadMats.Q(3:3:end,:)*dXdpstat;
    prevS.uxykp0 = zeros(size(prevS.uxynFxyn, 1), 2); %Tracks reference point of kp spring
    
    %% Modal Analysis
    [Vst, Wst] = eigs(dRstat, M, 10, 'SM');
    [Wst,si] = sort(sqrt(diag(Wst)));
    Vst = Vst(:,si)./sqrt(diag(Vst(:,si)'*M*Vst(:,si))');
    
    Vmi = Vst(:, mdi);
    Wmi = Wst(mdi);
    %% RQNM
    BB.Q = expdat.Q;
    BB.W = zeros(size(BB.Q));
    BB.Z = zeros(size(BB.Q));
    BB.D = zeros(size(BB.Q));
    BB.X = zeros(length(Xstat)+1, length(BB.Q));
    BB.R = zeros(length(Xstat)+1, length(BB.Q));
    
    %I am not worrying about derivatives with parameters too much right
    %now.
%     dWdp = zeros(1, numel(pars));
%     dZdp = zeros(1, numel(pars));
%     dDdp = zeros(1, numel(pars));
    
    ErrsWZD = zeros(3,1);
    dErrsdp = zeros(3, numel(pars));
    
    % LGL Quadrature Points
    [xi, wi] = LGLWT(Nqp, 0, 1);
    qul = zeros(size(xi));  aul = zeros(size(xi));
    qrl = zeros(size(xi));  arl = zeros(size(xi));
    
    
    for i=1:length(BB.Q)
        
        dDdp = zeros(1,numel(pars));   
        
        % Backbone Loading
        prev = prevS;
        
        %Backbone - jump straight to the first reversal point. 
        qbb = BB.Q(i);
        X0 = [Xstat+Vmi*qbb; Wmi^2];

        Xb = fsolve(@(X) WJRES_RQNM_LOOP([X; qbb], pars, prev, K, M, Xstat, dXdpstat, Fv, L, QuadMats.T, QuadMats.Q, CFUN, Npatches), X0, opt);
        [~, dRdX, ~, dRdp] = WJRES_RQNM_LOOP([Xb; qbb], pars, prev, K, M, Xstat, dXdpstat, Fv, L, QuadMats.T, QuadMats.Q, CFUN, Npatches);

        dXbdp = -dRdX\dRdp;

        %Update prev
        uxyn = [QuadMats.Q(1:3:end,:)*Xb(1:end-1), QuadMats.Q(2:3:end,:)*Xb(1:end-1), QuadMats.Q(3:3:end,:)*Xb(1:end-1)]; %WJ uxyn
        [Fxyn, dFxyn, dFxyndp, prev] = CFUN(uxyn, pars, prev);

        % Hysteretic Unloading
        prev.uxynFxyn = [uxyn, Fxyn];
        prev.dFxyn = dFxyn;
        prev.dFxyndp = dFxyndp;
        prev.duxdp = QuadMats.Q(1:3:end,:)*dXbdp(1:end-1, :);
        prev.duydp = QuadMats.Q(2:3:end,:)*dXbdp(1:end-1, :);
        prev.dundp = QuadMats.Q(3:3:end,:)*dXbdp(1:end-1, :);
%         prev.uxykp0 = prev.uxykp0; %included in the static prev updating if needed with penalty methods
    
        for k=1:(Nqp+1)
            qul(k) = BB.Q(i) - 2*BB.Q(i)*xi(k);
            X0 = [Xstat+Vmi*qul(k); Wmi^2];
            
            Xd = fsolve(@(X) WJRES_RQNM_LOOP([X; qul(k)], pars, prev, K, M, Xstat, dXdpstat, Fv, L, QuadMats.T, QuadMats.Q, CFUN, Npatches), X0, opt);
            [~, dRdX, ~, dRdp] = WJRES_RQNM_LOOP([Xd; qul(k)], pars, prev, K, M, Xstat, dXdpstat, Fv, L, QuadMats.T, QuadMats.Q, CFUN, Npatches);
                        
            %Update prev
%             uxyn = QuadMats.Q*reshape(L(1:MESH.Nn*3,:)*Xd(1:end-1), 3, [])';
            uxyn = [QuadMats.Q(1:3:end,:)*Xd(1:end-1), QuadMats.Q(2:3:end,:)*Xd(1:end-1), QuadMats.Q(3:3:end,:)*Xd(1:end-1)]; %WJ uxyn
            [Fxyn, dFxyn, dFxyndp, prev] = CFUN(uxyn, pars, prev);
            
            dXdp = -dRdX\dRdp;
            aul(k) = qul(k)*Xd(end);
            
            BB.D(i) = BB.D(i) - qul(k)*wi(k)*(2*BB.Q(i))*Xd(end);
            dDdp = dDdp - qul(k)*wi(k)*(2*BB.Q(i))*dXdp(end,:);
        end
        
        % Hysteretic Reloading
        prev.uxynFxyn = [uxyn, Fxyn];
        prev.dFxyn = dFxyn;
        prev.dFxyndp = dFxyndp;
        prev.duxdp = QuadMats.Q(1:3:end,:)*dXdp(1:end-1, :);
        prev.duydp = QuadMats.Q(2:3:end,:)*dXdp(1:end-1, :);
        prev.dundp = QuadMats.Q(3:3:end,:)*dXdp(1:end-1, :);
%         prev.uxykp0 = prev.uxykp0; %included in the static prev updating if needed with penalty methods
        
        
        for k=1:(Nqp+1)
            qrl(k) = -BB.Q(i) + 2*BB.Q(i)*xi(k);
            X0 = [Xstat+Vmi*qrl(k); Wmi^2];

            Xd = fsolve(@(X) WJRES_RQNM_LOOP([X; qrl(k)], pars, prev, K, M, Xstat, dXdpstat, Fv, L, QuadMats.T, QuadMats.Q, CFUN, Npatches), X0, opt);
            [R, dRdX, ~, dRdp] = WJRES_RQNM_LOOP([Xd; qrl(k)], pars, prev, K, M, Xstat, dXdpstat, Fv, L, QuadMats.T, QuadMats.Q, CFUN, Npatches);
                
            
            %Update prev
            uxyn = [QuadMats.Q(1:3:end,:)*Xd(1:end-1), QuadMats.Q(2:3:end,:)*Xd(1:end-1), QuadMats.Q(3:3:end,:)*Xd(1:end-1)]; %WJ uxyn
            [~, ~, ~, prev] = CFUN(uxyn, pars, prev);
            
            dXdp = -dRdX\dRdp;
            arl(k) = qrl(k)*Xd(end);
            
            BB.D(i) = BB.D(i) + Xd(end)*qrl(k)*wi(k)*(2*BB.Q(i));
            dDdp = dDdp + qrl(k)*wi(k)*(2*BB.Q(i))*dXdp(end,:);
        end
        
        % Processing
        BB.W(i) = sqrt(Xb(end));
        BB.Z(i) = BB.D(i)/(2*pi*(BB.Q(i)*BB.W(i))^2) + expdat.Z(1);
        
        dWdp = dXbdp(end,:)/(2*BB.W(i));
%         dZdp =
%         (dDdp-BB.D(i)/Xb(end)*dXbdp(end,:))/(2*pi*(BB.Q(i)*BB.W(i))^2);
%         %This was just a different way of writing the same thing.  
        dZdp = (dDdp - 2*BB.D(i)/BB.W(i)*dWdp)/(2*pi*(BB.Q(i)*BB.W(i))^2);
        
        if isreal(BB.W(i))
            ErrsWZD(1) = ErrsWZD(1) + ((BB.W(i)-expdat.W(i))/expdat.W(i))^2;
            dErrsdp(1,:) = dErrsdp(1,:) + 2*((BB.W(i)-expdat.W(i))/expdat.W(i)^2)*dWdp;
        else
            ErrsWZD(1) = ErrsWZD(1) + 1;
            dErrsdp(1,:) = dErrsdp(1,:) + 0;
        end
        
        if BB.Z(i)>0
            ErrsWZD(2) = ErrsWZD(2) + ((log(BB.Z(i))-log(expdat.Z(i)))/log(expdat.Z(i)))^2;
            ErrsWZD(3) = ErrsWZD(3) + ((log(BB.D(i))-log(expdat.D(i)))/log(expdat.D(i)))^2;
            dErrsdp(2,:) = dErrsdp(2,:) + 2/BB.Z(i)*((log(BB.Z(i))-log(expdat.Z(i)))/log(expdat.Z(i))^2)*dZdp;
            dErrsdp(3,:) = dErrsdp(3,:) + 2/BB.Z(i)*((log(BB.D(i))-log(expdat.D(i)))/log(expdat.D(i))^2)*dDdp;
        else
            ErrsWZD(2) = ErrsWZD(2) + 1;
            ErrsWZD(3) = ErrsWZD(3) + 1;
            dErrsdp(2,:) = dErrsdp(2,:) + 0;
            dErrsdp(3,:) = dErrsdp(3,:) + 0;
        end
        
        BB.X(:, i) = Xd;
        BB.R(:, i) = R;
        fprintf('%d...', i);
%         toc
    end
%     fprintf('\nCompleted RQNM....................................');
    
    ErrsWZD = ErrsWZD/length(BB.Q);
    dErrsdp = dErrsdp/length(BB.Q);
    
    ErrsWZD = ErrsWZD([1 2], :);
    dErrsdp = dErrsdp([1 2], :);
    
    % Log-Scale Errors
    dErrsdp = (dErrsdp./ErrsWZD)/log(10);
    ErrsWZD = log10(ErrsWZD);
    
    if ~isreal(ErrsWZD)
        disp('complex');
    end    
end