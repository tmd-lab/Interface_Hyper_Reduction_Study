function [F,J] = ZTECONTACTFORCEJAC(contactfunc,U,MESH,p0,pars,QuadMats,prev, varargin)
%ZTE_QP_CONTACTFORCEJAC Returns the force and Jacobian after
%applying the traction-based contact model.
% USAGE:
%	[F,J] = ZTE_QP_CONTACTFORCEJAC(contactfunc,U,MESH,p0,pars,QuadMats,prev);
% INPUTS:
%   contactfunc	: Function handle
%   U		: Nodal Displacements
%   MESH	: Mesh structure 
%   p0		: Static nodal pressures 
%   pars	: friction parameters 
%   QuadMats	: Quadrature Q & T matrices 
%   prev	: Previous states at nodal/quadrature points
%   		  (interpolated to quadrature locations if nodal)
% OUTPUTS:
%   F		: Nodal Force vector
%   J		: Nodal Jacobian matrix

    Ndof 	= length(U);
    Nq      = size(QuadMats.Q,1);
    uxyn_qp	= QuadMats.Q*reshape(U(1:(MESH.Nn*3)), 3, MESH.Nn)'; % 3 dofs at
                                                  % quadrature
                                                  % locations

    % Interpolating previous and static states expressed at nodal
    % points to the quadrature locations.
    if size(prev.uxyntxyn,1)==MESH.Nn
        prev.uxyntxyn	= QuadMats.Q*prev.uxyntxyn;
    end
    if length(p0)==MESH.Nn
        p0 	= QuadMats.Q*p0;
    end
    if length(p0)==1
        p0 = ones(size(QuadMats.Q, 1),1)*p0;
    end
    [Ptx,Pty,Pn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun] = ...
        contactfunc(uxyn_qp, p0, pars, prev.uxyntxyn);
    
    % Integrating & Assembling - Forces
    F = zeros(Ndof,1);
    F(1:(3*MESH.Nn)) = reshape((QuadMats.T*[Ptx Pty Pn])',3*MESH.Nn,1);
    
    % Integrating & Assembling - Jacobians
    J = sparse(Ndof,Ndof);
    
    J(1:3:(MESH.Nn*3), 1:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtxdux,1,size(QuadMats.Q,2)).*QuadMats.Q);
    J(1:3:(MESH.Nn*3), 2:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtxduy,1,size(QuadMats.Q,2)).*QuadMats.Q);
    J(1:3:(MESH.Nn*3), 3:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtxdun,1,size(QuadMats.Q,2)).*QuadMats.Q);
    
    J(2:3:(MESH.Nn*3), 1:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtydux,1,size(QuadMats.Q,2)).*QuadMats.Q);
    J(2:3:(MESH.Nn*3), 2:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtyduy,1,size(QuadMats.Q,2)).*QuadMats.Q);
    J(2:3:(MESH.Nn*3), 3:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtydun,1,size(QuadMats.Q,2)).*QuadMats.Q);
    
    J(3:3:(MESH.Nn*3), 3:3:(MESH.Nn*3)) = QuadMats.T*(repmat(dtndun,1,size(QuadMats.Q,2)).*QuadMats.Q);
    
    if nargin==8
        F = varargin{1}*F;
        J = varargin{1}*J;
    end
% %     Checking for pinning
%     if isfield(MESH, 'Pstiff')
%         for b=1:length(MESH.Ndbs)
%             Ndb = length(MESH.Ndbs{b});
%             xy = U((MESH.Ndbs{b}-1)*3+(1:2));
%             
%             xyC = mean(xy);
%             rC = vecnorm(xyC);
%             xyhatC = xyC/rC;
%             
%             drCdx = xyhatC(1)*ones(1,Ndb)/Ndb;
%             drCdy = xyhatC(2)*ones(1,Ndb)/Ndb;
%             
%             if rC>MESH.clrn
%                 f = MESH.Pstiff*(rC-MESH.clrn);
%                 fxy = f*xyhatC;
%             
%                 F((MESH.Ndbs{b}-1)*3+1) = F((MESH.Ndbs{b}-1)*3+1) + fxy(1)/Ndb;
%                 F((MESH.Ndbs{b}-1)*3+2) = F((MESH.Ndbs{b}-1)*3+2) + fxy(2)/Ndb;
%             
%                 J((MESH.Ndbs{b}-1)*3+1, (MESH.Ndbs{b}-1)*3+1) = J((MESH.Ndbs{b}-1)*3+1, (MESH.Ndbs{b}-1)*3+1) + ...
%                     MESH.Pstiff*(1-MESH.clrn/rC*(1-xyhatC(1)^2))*ones(Ndb)/Ndb^2;
%                 J((MESH.Ndbs{b}-1)*3+1, (MESH.Ndbs{b}-1)*3+2) = J((MESH.Ndbs{b}-1)*3+1, (MESH.Ndbs{b}-1)*3+2) + ...
%                     MESH.Pstiff*(MESH.clrn/rC*xyhatC(1)*xyhatC(2))*ones(Ndb)/Ndb^2;
%                 
%                 J((MESH.Ndbs{b}-1)*3+2, (MESH.Ndbs{b}-1)*3+1) = J((MESH.Ndbs{b}-1)*3+2, (MESH.Ndbs{b}-1)*3+1) + ...
%                     MESH.Pstiff*(MESH.clrn/rC*xyhatC(1)*xyhatC(2))*ones(Ndb)/Ndb^2;                
%                 J((MESH.Ndbs{b}-1)*3+2, (MESH.Ndbs{b}-1)*3+2) = J((MESH.Ndbs{b}-1)*3+2, (MESH.Ndbs{b}-1)*3+2) + ...
%                     MESH.Pstiff*(1-MESH.clrn/rC*(1-xyhatC(2)^2))*ones(Ndb)/Ndb^2;
%             end
%         end
%     end
end