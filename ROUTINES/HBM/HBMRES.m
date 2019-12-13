function [R,dRdX,dRdw,dRdXw] = HBMRES(X,M,C,K,Fl,L,h,QuadMats,contactfun, varargin)
%HBMRES Returns the residual and derivatives for the harmonic
%balance system
% USAGE:
%	[R,dRdXa] = HBMRES(X,M,C,K,Fl,L,h,MESH,QuadMats,contactfun);
% INPUTS:
%   X,MESH,QuadMats,No 
% OUTPUTS:
%   R,dRdXa

    if nargin>=11
        X = varargin{2}.*X;
    end

    [E,dEdw] 	= HARMONICSTIFFNESS(M,C,K,X(end),h);
    Nd		= size(M,1);                % Model DOF's
    Nhc 	= length(X(1:end-1))/Nd;    % No of harmonic components
    Nc          = size(QuadMats.Qst, 2);    % No of nonlinear components
    Nrest	= size(L,1)-Nc;             % No of generalized components
    Nph     = Nc+Nrest;               % Number of physical DOF's
    Nq      = size(QuadMats.Q,1);       % Number of quadrature points
    
    Xds     = L*reshape(X(1:end-1),Nd,Nhc);         % Physical DOF's
    Xqp     = (QuadMats.Qst*Xds(1:Nc,:))';   % Quadrature points
                                                    % [qp1 qp2 qp3...]
    
    [Tx,Ty,Tn,Jxx,Jxy,Jxn,Jyx,Jyy,Jyn,Jnn] = ...
        contactfun(Xqp(:,1:Nq),Xqp(:,Nq+(1:Nq)),Xqp(:,2*Nq+(1:Nq)), h);
    
    tmp = kron(speye(Nq), ones(Nhc));
    ltmp = logical(tmp);
    tmp(ltmp) = Jxx; Jxx = tmp;
    tmp(ltmp) = Jxy; Jxy = tmp;
    tmp(ltmp) = Jxn; Jxn = tmp;
    
    tmp(ltmp) = Jyx; Jyx = tmp;
    tmp(ltmp) = Jyy; Jyy = tmp;
    tmp(ltmp) = Jyn; Jyn = tmp;
    
    tmp(ltmp) = Jnn; Jnn = tmp;
%     clear tmp ltmp

    Jaux = sparse([Jxx Jxy Jxn;
                   Jyx Jyy Jyn;
                   sparse(Nq*Nhc,2*Nq*Nhc) Jnn]);
%     clear Jxx Jxy Jxn Jyx Jyy Jyn Jnn

    rcs = 1:size(Jaux,1);
    for i=1:Nhc
        rcs((i-1)*Nq*3+(1:Nq*3)) = i:Nhc:(Nq*3*Nhc);
    end
    Jro = Jaux(rcs,rcs);
%     clear Jaux

    Fnl = reshape(L'*[QuadMats.Tst*[Tx Ty Tn]'; zeros(Nrest,Nhc)],Nd*Nhc,1);
    
    Jro = kron(speye(Nhc),QuadMats.Tst)*Jro*kron(speye(Nhc),QuadMats.Qst);
    Jnlp = sparse(Nph*Nhc,Nph*Nhc);
%     % PREV
%     rcs = zeros((MESH.Nn*3)*Nhc,1);
%     for i=1:Nhc
%         rcs((i-1)*(MESH.Nn*3)+(1:(MESH.Nn*3))) = (i-1)*(Nph)+(1:MESH.Nn*3);
%     end
    % CUR
    rcs = zeros(Nc*Nhc,1);
    for i=1:Nhc
        rcs((i-1)*Nc+(1:Nc)) = (i-1)*(Nph)+(1:Nc);
    end
%     % Vectorized Form
%     rcs = kron(((1:Nhc)-1)*Nph, ones(1,MESH.Nn*3)) +
%     kron(ones(1,Nhc),(1:MESH.Nn*3));
    Jnlp(rcs,rcs) = Jro;
    tmp = sparse(kron(speye(Nhc),L));
    %%%%%%%%% TRADE OFF HERE FOR SPEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tic
%     Jnl = tmp'*Jnlp*tmp;
%     toc
%     tic
    Jnl = zeros(Nd*Nhc,Nd*Nhc);
    for i=1:Nhc
        for j=1:Nhc
            Jnl((i-1)*Nd+(1:Nd),(j-1)*Nd+(1:Nd)) = L'*Jnlp((i-1)*Nph+(1:Nph),(j-1)*Nph+(1:Nph))*L;
        end
    end
%     toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear Tx Ty Tn Jro
    
    R		= E*X(1:end-1) + Fnl - Fl;
    dRdX    = E+Jnl;
    dRdw    = dEdw*X(1:end-1);

    if nargin>=11
        dRdX = dRdX.*varargin{2}(1:end-1);
        dRdw = dRdw*varargin{2}(end);
    end
    if nargin>=10
        R = varargin{1}.*R;
        dRdX = varargin{1}.*dRdX;
        dRdw = varargin{1}.*dRdw;
    end
    
    dRdXw 	= [dRdX dRdw];
end
