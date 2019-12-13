function [R,dRdX,dRdw,dRdXa] = EPMCRES(X,M,C,K,Fl,L,h,QuadMats,contactfun, varargin)
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
    
    Nd		= size(M,1);                % Model DOF's
    Nhc 	= sum(h==0) + 2*sum(h~=0);    % No of harmonic components
    Nc      = size(QuadMats.Qst, 2);    % No of nonlinear components
    Nrest	= size(L,1)-Nc;             % No of generalized components
    Nph     = Nc+Nrest;               % Number of physical DOF's
    Nq      = size(QuadMats.Q,1);       % Number of quadrature points
    
    [E,dEdw] 	= HARMONICSTIFFNESS(M,C-X(end-1)*M,K,X(end-2),h);
    dEdxi       = HARMONICSTIFFNESS(zeros(size(M)), -M, zeros(size(K)), X(end-2), h);
    
    amp     = 10.^X(end);
    dampdXe = log(10)*amp;
    Asc     = kron([1; amp*ones(Nhc-1,1)], ones(Nd,1));
    dXdamp  = kron([0; ones(Nhc-1,1)], ones(Nd,1));
    
	Xscl    = Asc.*X(1:end-3);
    Xqp     = (QuadMats.Qst*(L(1:Nc,:)*reshape(Xscl,Nd,Nhc)))';   % Quadrature points
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

    Jaux = sparse([Jxx Jxy Jxn;
                   Jyx Jyy Jyn;
                   sparse(Nq*Nhc,2*Nq*Nhc) Jnn]);

    rcs = 1:size(Jaux,1);
    for i=1:Nhc
        rcs((i-1)*Nq*3+(1:Nq*3)) = i:Nhc:(Nq*3*Nhc);
    end
    Jro = Jaux(rcs,rcs);

    Fnl = reshape(L'*[QuadMats.Tst*[Tx Ty Tn]'; zeros(Nrest,Nhc)],Nd*Nhc,1);

    Jro = kron(speye(Nhc),QuadMats.Tst)*Jro*kron(speye(Nhc),QuadMats.Qst);
    Jnlp = sparse(Nph*Nhc,Nph*Nhc);

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

    R		= E*Xscl + Fnl - [Fl(1:Nd); 0*Fl(Nd+1:end)];
    dRdX    = (E+Jnl).*Asc;
    dRdw    = dEdw*Xscl;
    dRdxi   = dEdxi*Xscl;
    
    % Phase Condition for EPMC
    PFl = Fl;   PFl(1:Nd) = 0;
    R   = [R; ...
           X(Nd+(1:Nd))'*M*X(Nd+(1:Nd))+X(2*Nd+(1:Nd))'*M*X(2*Nd+(1:Nd))-1.0; ...
           PFl'*X(1:end-3)];

    dRdX = [dRdX dRdw dRdxi; ...
            zeros(1,Nd) 2*X(Nd+(1:Nd))'*M 2*X(2*Nd+(1:Nd))'*M zeros(1,(Nhc-3)*Nd) 0 0; ...
            PFl' 0 0];
    dRda = [(E+Jnl)*(dXdamp.*X(1:end-3))*dampdXe; ...
            0; 0];
    
    % Conditioning
    if nargin>=10  % Left Pre Conditioner
        R = varargin{1}.*R;
        dRdX = varargin{1}.*dRdX;
        dRda = varargin{1}.*dRda;
    end
    if nargin>=11
        dRdX = dRdX.*varargin{2}(1:end-1);
        dRda = dRda*varargin{2}(end);
    end
    
    dRdXa 	= [dRdX dRda];
end
