function [Mr,Kr,TFM] = HCBREDUCE(M,K,bdofs,ncomp, varargin)
%HCBREDUCE Returns the Fixed Interface Hurty/Craig-Bampton Reduced
%Mass and Stiffness Matrices for a given symmetric system.
% USAGE:
%	[Mr,Kr,TFM] = HCBREDUCE(M,K,bdofs,ncomp);
% INPUTS:
%   M,K		: NdofxNdof mass & stiffness matrices 
%   bdofs	: Nbx1 set of boundary/retained DOF's
%   ncomp	: 1x1 Number of fixed boundary modes
%   varargin{1} : scalar float. Use null space transform in HCB solution to
%                   eliminate any singular values with S(i)/S(1) < varargin{1}
%               If not included in input, then no null space transform 
%               is done for the matrix solution in HCB.
%               Recommended value based on preliminary tests would be 1e-10
% OUTPUTS:
%   Mr,Kr	: (Ndof-Nb)x(Ndof-Nb) Reduced mass & stiffness
%   		  matrices
%   TFM		: Ndofx(Ndof-Nb) Transformation matrix

    Ndof = size(M,1);
    bdofs = reshape(bdofs, 1, length(bdofs));
    idofs = setdiff(1:Ndof,bdofs);

    Kbb = 0.5*(K(bdofs,bdofs)+K(bdofs,bdofs)');
    Kbi = 0.5*(K(bdofs,idofs)+K(idofs,bdofs)');
    Kii = 0.5*(K(idofs,idofs)+K(idofs,idofs)');

    Mbb = 0.5*(M(bdofs,bdofs)+M(bdofs,bdofs)');
    Mbi = 0.5*(M(bdofs,idofs)+M(idofs,bdofs)');
    Mii = 0.5*(M(idofs,idofs)+M(idofs,idofs)');
    clear K M

    if nargin==4  % No expected nulls
        PhiC = -Kii\(Kbi'); % Constraint modes
    else
%         [Vii, Dii] = eigs(Kii, varargin{1}*3, 'SM');  [Dii,si] = sort(diag(Dii));
%         Vii = Vii(:, si);
%         L = null(Vii(:, 1:varargin{1})');

        % Previous algorithm
%         L = null( null(full(Kii))' );

        % Manual null space in one SVD call
        [U, S, V] = svd(full(Kii)); % Kii = U * S * V'
        num_drop = sum(diag(S) / S(1,1) < varargin{1});
        
        L = V(:, 1:end-num_drop); % Columns of V are a basis for X in Kii*X

        PhiC = -L* ((L'*Kii*L)\(L'*Kbi'));
    end
    try
        [PhiN,Dfb] = eigs(Kii,Mii,ncomp*2,'SM');
    catch me
        fprintf('Attemping full eigen solution\n');
        [PhiN,Dfb] = eig(full(Kii),full(Mii));
    end
    [Dfb,si] = sort(diag(Dfb));
    PhiN = PhiN(:,si);  PhiN = PhiN./sqrt(diag(PhiN'*Mii*PhiN)');
    Dfb = Dfb((1:ncomp));
    PhiN = PhiN(:,(1:ncomp));
    
    Mcc = Mbb+Mbi*PhiC+PhiC'*Mbi'+PhiC'*Mii*PhiC;
    Mcn = (Mbi+PhiC'*Mii)*PhiN;
  
    Mr = [0.5*(Mcc+Mcc') Mcn; Mcn' eye(ncomp)];%PhiN'*Mii*PhiN];
    Kr = [Kbb+0.5*(Kbi*PhiC+PhiC'*Kbi') zeros(size(Mcn)); 
        zeros(size(Mcn')) diag(Dfb)];%PhiN'*Kii*PhiN];
    TFM = sparse([eye(length(bdofs)) zeros(length(bdofs),ncomp);
           PhiC PhiN]);
    TFM([bdofs idofs], :) = TFM;
end