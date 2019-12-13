function [F, J] = GENELCONTACTFORCEJAC(contactfunc,U,MESH,p0,pars,QuadMats,prev)
	Ndof 	= length(U);
    Nqdof   = size(QuadMats.Q,1);
    uxyn_qp = reshape(QuadMats.Q*U(1:MESH.Nn*MESH.dpn), 3, Nqdof/3)';

    if length(p0)==1
        p0 = ones(Nqdof/3,1)*p0;
    end
    [Ptx,Pty,Pn,dtxdux,dtxduy,dtxdun,dtydux,dtyduy,dtydun,dtndun] = ...
        contactfunc(uxyn_qp, p0, pars, prev.uxyntxyn);
    
    % Integrating & Assembling - Forces
    F = zeros(Ndof,1);
    F(1:(MESH.Nn*MESH.dpn)) = QuadMats.T*reshape([Ptx Pty Pn]', Nqdof, 1);
    % Integrating & Assembling - Jacobians
    J = sparse(Ndof,Ndof);
    
    Jtmp = zeros(Nqdof, Nqdof);
    Jtmp(1:3:end, 1:3:end) = diag(dtxdux);
    Jtmp(1:3:end, 2:3:end) = diag(dtxduy);
    Jtmp(1:3:end, 3:3:end) = diag(dtxdun);
    
    Jtmp(2:3:end, 1:3:end) = diag(dtydux);
    Jtmp(2:3:end, 2:3:end) = diag(dtyduy);
    Jtmp(2:3:end, 3:3:end) = diag(dtydun);
    
    Jtmp(3:3:end, 3:3:end) = diag(dtndun);
    
    J(1:MESH.Nn*MESH.dpn, 1:MESH.Nn*MESH.dpn) = QuadMats.T*Jtmp*QuadMats.Q;
end