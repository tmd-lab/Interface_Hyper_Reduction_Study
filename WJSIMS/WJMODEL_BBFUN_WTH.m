function [ErrsWZD, dErrsdp, BB] = WJMODEL_BBFUN_WTH(pars, mdi, expdat, K, M, R, X0, Fv, L, Txyntxyn, Qxyntxyn, CFUN, Npatches, Nqp, opt)



    % Prestress
    [Xstat, eflg] = fsolve(@(X) WJRES_STR_WTH([X; 0.0], pars, K, Fv*0, Fv, L, Txyntxyn, Qxyntxyn, CFUN, Npatches), X0, opt);
    if eflg<0
        disp(pars)
        error('Non-convergent!');
    end
    [~,dRstat,~,dRdpstat] = WJRES_STR_WTH([Xstat; 0.0], pars, K, Fv*0, Fv, L, Txyntxyn, Qxyntxyn, CFUN, Npatches);
    dXdpstat = -dRstat\dRdpstat;
	
%     prev.uxyn = [Qxyn(1:3:end,:)*Xstat, Qxyn(2:3:end,:)*Xstat, Qxyn(3:3:end,:)*Xstat*0];
%     [prev.Fxyn, dFxynstat, prev.dFxyndpar] = CFUN(prev.uxyn, pars);
    
    prev.uxyntxyn = [Qxyntxyn(1:6:end,:)*Xstat, Qxyntxyn(2:6:end,:)*Xstat, Qxyntxyn(3:6:end,:)*Xstat Qxyntxyn(4:6:end,:)*Xstat, Qxyntxyn(5:6:end,:)*Xstat, Qxyntxyn(6:6:end,:)*Xstat];
    [prev.Fxyntxyn, dFxyntxynstat, prev.dFxyntxyndpar] = CFUN(prev.uxyntxyn, pars);
    prev.uxyntxyn(:, 3) = prev.uxyntxyn(:, 3)*0;
    prev.Fxyntxyn(:, 3) = prev.Fxyntxyn(:, 3)*0;
    
	if size(dFxyntxynstat,3)==1
       tmp = dFxyntxynstat;
       dFxyntxynstat = zeros(size(dFxyntxynstat,1),size(dFxyntxynstat,2),3);
       dFxyntxynstat(:,1,1) = tmp(:,1);
       dFxyntxynstat(:,2,2) = tmp(:,2);
       dFxyntxynstat(:,3,3) = tmp(:,3); % Np x uxyn x fxyn
    end
	dFxyntxynstat = num2cell(permute(dFxyntxynstat,[3,2,1]), [2,1]);
    prev.dFxyntxyndXdpar = blkdiag(dFxyntxynstat{:})*Qxyntxyn*dXdpstat;
    
%     prev.dFxyndXdpar = reshape(dFxynstat', Npatches*3, 1).*Qxyn*dXdpstat;
    
%     prev.uxyn = prev.uxyn*0;
%     prev.Fxyn = prev.Fxyn*0;
%     prev.dFxyndpar = prev.dFxyndpar*0;
%     prev.dFxyndXdpar = prev.dFxyndXdpar*0;
    %% Modal Analysis
    [Vst, Wst] = eigs(dRstat, M, 10, 'SM');
    Wst = sqrt(diag(Wst));
    Vst = Vst./sqrt(diag(Vst'*M*Vst)');

    Vm = Vst(:, mdi);
    Wm = Wst(mdi);
    %% QSNMA
    BB.Q = expdat.Q;
    BB.W = zeros(size(BB.Q));
    BB.Z = zeros(size(BB.Q));
    BB.D = zeros(size(BB.Q));
    BB.Phi = zeros(length(Vm), length(BB.Q));
    
    dWdp = zeros(1, length(pars));
    dDdp = zeros(1, length(pars));
    dZdp = zeros(1, length(pars));
    
    ErrsWZD = zeros(3, 1);
    dErrsdp = zeros(3, length(pars));
    
    % LGL Quadrature points
    [xi, wi] = LGLWT(Nqp, 0, 1);
    qr = zeros(size(xi));
    lr = qr;
    
    % Initialization for accumulators
    Q0 = 0;
    Di = 0;
    dDdpi = zeros(1, length(pars));
    Vmi = Vm;
    Wmi = Wm;
    
    rmp = 0;  % Number of points to remove
    for i=1:length(BB.Q)
        if i>1
            Q0 = BB.Q(i-1);
        end
        for k=2:(Nqp+1)  % Evaluate Nqp points for quadrature
            qr(k) = Q0 + (BB.Q(i)-Q0)*xi(k);
            
            if qr(k)==0
                lr(k) = Wmi^2;
                continue;
            end
            
            X0 = [Xstat+Vmi*qr(k); Wmi^2];
            Xd = fsolve(@(X) WJRES_RQNM_WTH([X; qr(k)], pars, prev, K, M, Xstat, dXdpstat, ...
                Fv, L, Txyntxyn, Qxyntxyn, CFUN, Npatches), X0, opt);
            lr(k) = Xd(end);
            [~, dRdX, ~, dRdp] = WJRES_RQNM_WTH([Xd; qr(k)], pars, prev, K, M, Xstat, dXdpstat, Fv, ...
                L, Txyntxyn, Qxyntxyn, CFUN, Npatches);
            dXdp = -dRdX\dRdp;
            
            % LGL Quadrature
            Di = Di + (BB.Q(i)-Q0)*wi(k)*qr(k)*abs(Xd(end));
            dDdpi = dDdpi + (BB.Q(i)-Q0)*wi(k)*qr(k)*dXdp(end,:);
        end
        Vmi = (Xd(1:end-1)-Xstat)/BB.Q(i);
        Wmi = sqrt(Xd(end));
        
        BB.D(i) = -4*BB.Q(i)^2*Xd(end) + 8*Di;
        dDdp = -4*BB.Q(i)^2*dXdp(end,:) + 8*dDdpi;
        
        BB.W(i) = sqrt(Xd(end));
        BB.Z(i) = BB.D(i)/(2*pi*(BB.Q(i)*BB.W(i))^2) + expdat.Z(1)*0;
        BB.Phi(:, i) = (Xd(1:end-1)-Xstat)/BB.Q(i);
        
        dWdp = dXdp(end,:)/(2*BB.W(i));
        dZdp = (dDdp - 2*BB.D(i)/BB.W(i)*dWdp)/(2*pi*(BB.Q(i)*BB.W(i))^2);
        
        if isreal(BB.W(i))
            ErrsWZD(1) = ErrsWZD(1) + ((BB.W(i)-expdat.W(i))/expdat.W(i))^2;
            dErrsdp(1,:) = dErrsdp(1,:) + 2*((BB.W(i)-expdat.W(i))/expdat.W(i)^2)*dWdp;
        else
            ErrsWZD(1) = ErrsWZD(1) + 1;
            dErrsdp(1,:) = dErrsdp(1,:) + 0;
        end
        
%         ErrsWZD(2) = ErrsWZD(2) + ((BB.Z(i)-expdat.Z(i))/expdat.Z(i))^2;
%         ErrsWZD(3) = ErrsWZD(3) + ((BB.D(i)-expdat.D(i))/expdat.D(i))^2;
%         dErrsdp(2,:) = dErrsdp(2,:) + 2*((BB.Z(i)-expdat.Z(i))/expdat.Z(i)^2)*dZdp;
%         dErrsdp(3,:) = dErrsdp(3,:) + 2*((BB.D(i)-expdat.D(i))/expdat.D(i)^2)*dDdp;
        
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
        fprintf('%d ',i);
    end
    fprintf('\n');
    ErrsWZD = ErrsWZD/length(BB.Q);
    dErrsdp = dErrsdp/length(BB.Q);
    
    ErrsWZD = ErrsWZD([1 2],:);
    dErrsdp = dErrsdp([1 2],:);

    % Log-Scale Errors
    dErrsdp = (dErrsdp./ErrsWZD)/log(10);
    ErrsWZD = log10(ErrsWZD);
    
    if ~isreal(ErrsWZD)
        disp('complex');
    end
end
