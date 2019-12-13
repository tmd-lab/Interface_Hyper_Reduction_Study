function [Dtot,dx,dy,dn,Dx,Dy,Dn] = ZTEDISSIPATION_GEN(U,P,QuadMats,MESH)
%ZTEDISSIPATION Calculates the energy dissipated over a cycle of
%data provided
% USAGE:
%	[Dtot,D] = ZTEDISSIPATION(U,T,MESH);
% INPUTS:
%   U		: (NnxNcyc) Displacements 
%   T		: (NnxNcyc) Tractions 
%   MESH        : Structure containing interface mesh information
%       Nn      : Number of interface nodes
%       Nds     : Nnx2 [x y] coordinates of interface nodes
%       Ne      : Number of interface elements
%       Ne_Tri	: Number of interface triangular elements
%       Ne_Quad	: Number of interface quadrilateral elements
%       Tri     : Ne_Trix4 [eid n1 n2 n3]
%       Quad	: Ne_Quadx4 [eid n1 n2 n3 n4]
%  No		:     
% OUTPUTS:
%   Dtot	: (1x1) Total Integrated dissipation 
%   D		: (Nnx1) Dissipation distribution estimate (inaccurate)

    Nn = MESH.Nn;
    Nq = size(QuadMats.Q,1)/3;
    
    Uxyn = QuadMats.Q*U(1:(MESH.dpn*MESH.Nn), :);
    Ux = Uxyn(1:3:end, :);
    Uy = Uxyn(2:3:end, :);
    Un = Uxyn(3:3:end, :);
    
    Px = P(1:3:(3*Nq),:);
    Py = P(2:3:(3*Nq),:);
    Pn = P(3:3:(3*Nq),:);

    % Dissipation fluxes
    dx = dissfunc(Ux,Px);
    dy = dissfunc(Uy,Py);
    dn = dissfunc(Un,Pn);
    
    % Integrated nodal dissipation weights
    Dx = QuadMats.T(:,1:3:end)*dx;
    Dy = QuadMats.T(:,2:3:end)*dy;
    Dn = QuadMats.T(:,3:3:end)*dn;
    
    Dtot 	= sum(Dx+Dy+Dn*0);  % Don't consider normal direction dissipation
end

function [diss] = dissfunc(U,T)
    Np = size(U,1);
    % Trapz
    % for i=1:Np
    %    diss(i) = trapz(U(i,[0:end 1]),T(i,[1:end 1]));
    % end
    % Gauss/Shoelace Formula
%     for i=1:Np
%         diss(i) = 0.5*sum(U(i,:).*T(i,[2:end 1]) - U(i,[2:end 1]).*T(i,:));
%     end
    % Clockwise Shoelace
    diss = 0.5*(sum(U(:,[2:end 1]).*T, 2) - sum(U.*T(:,[2:end 1]), 2));
end
