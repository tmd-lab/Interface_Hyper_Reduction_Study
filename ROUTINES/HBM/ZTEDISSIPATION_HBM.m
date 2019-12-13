function [Dtot,Dx,Dy,Dn,dx,dy,dn] = ZTEDISSIPATION_HBM(X,L,MESH,QuadMats,h,Nhc,contactfun)
%ZTEDISSIPATION_HBM Returns the dissipation estimate for different
%forcing amplitudes given the harmonic balance solution
% USAGE:
%	[Dtot,Dx,Dy,Dn,dx,dy,dn] = ZTEDISSIPATION_HBM(X,L,MESH,QuadMats,h,Nhc,contactfun);
% INPUTS:
%   X,L,MESH,QuadMats,h,Nhc,contactfun 
% OUTPUTS:
%   Dtot,Dx,Dy,Dn,dx,dy,dn

    Np	= size(X,2);
    Nph	= size(L,1);
    Nq	= size(QuadMats.Q,1);
    % Quadrature Points
    Xqp = permute(reshape( kron(eye(Nhc), ...
                                QuadMats.Qst*L(1:MESH.Nn*3, :))*...
                           X(1:end-1, :), Nq*3, Nhc,Np), [2 1 3]);
    Xqx	= Xqp(:,1:Nq,:);
    Xqy = Xqp(:,Nq+(1:Nq),:);
    Xqn = Xqp(:,2*Nq+(1:Nq),:);
    clear Xqp
    
    Tx = zeros(size(Xqx));
    Ty = zeros(size(Xqx));
    Tn = zeros(size(Xqx));    
    
    Dm	= HARMONICSTIFFNESS(0,1,0,1.0,h);
    
    dx = zeros(Nq,Np);
    dy = dx;
    dn = dx;
    parfor i=1:Np
        [Tx(:,:,i),Ty(:,:,i),Tn(:,:,i)]	= ...
            contactfun(Xqx(:,:,i), Xqy(:,:,i), Xqn(:,:,i), h);
        dx(:,i)		= sum( Dm*Xqx(:,:,i).*Tx(:,:,i)*pi );
        dy(:,i)		= sum( Dm*Xqy(:,:,i).*Ty(:,:,i)*pi );
        dn(:,i)		= sum( Dm*Xqn(:,:,i).*Tn(:,:,i)*pi );
        disp([num2str(i) '/' num2str(Np) ' DONE.'])
    end
    clear Tx Ty Tn
    Dx = QuadMats.T*dx;
    Dy = QuadMats.T*dy;
    Dn = QuadMats.T*dn;
    Dtot = sum(Dx+Dy+Dn);
end