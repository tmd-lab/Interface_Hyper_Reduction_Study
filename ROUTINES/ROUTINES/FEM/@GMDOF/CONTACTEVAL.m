function [Fn, Z, dFndUn, dFndUnd, dFndPars, Fs_qp] = CONTACTEVAL(m, Un, Z, Und, Pars, varargin)
%CONTACTEVAL evaluates the given hysteretic contact forces at quadrature points and integrates them to provide nodal forces
%
% USAGE:
% ------
%   MESH.SETCFUN(fcont);  % Intialize
%   [Fn, Z, dFndUn, dFndUnd, dFndPars] = MESH.CONTACTEVAL(Un, Und, Z, Pars)
% INPUTS:
% -------
%   (Initialization)  
%   fcont	: Contact function handle of the form,
%		 	[fxyn, z, DfxynDuxyn, DfxynDuxynd, DfxynDpars] = fcont(uxyn, z, uxynd, Pars);
%   (Inputs)  
%   Un		: (MESH.Nn,1) Nodal displacements
%   Z		: (Nz,MESH2D.Ne*MESH2D.Nq^2)
%   Und		: (MESH.Nn,1) Nodal velocities  
%   Pars	: (Npars,1)  Parameters array
%   pA(optional): (xxx,Npars) Parameter access array. fcont will be called using
%			  reshape(pA*Pars, [], MESH.Ne*MESH.Nq^2) as Pars
%   pU(optional): (xxx,length(Un)) Projection/Mask for Un such that the Un that will be
%			  used for the calculations is mU*Un.
% OUTPUTS:
% --------
%   Fn		:
%   Z		:
%   dFndUn	:
%   dFndUnd	:
%   dFndPars	:

  if isempty(m.fcont)
    error('Contact force function not set');
  end

  if length(varargin)>=1
    pA = varargin{1};
    Parsc = reshape(pA*Pars, [], m.Nn);
  else
    pA = speye(length(Pars));
    Parsc = reshape(Pars, [], m.Nn);
  end
  Npars = size(Parsc, 1);

  if length(varargin)>=2
    pU = varargin{2};
    Un = pU*Un;
    Und = pU*Und;
  end

  Nu = length(Un);     % Number of dofs in U vector

  Us_qp = reshape(m.Qm*Un, m.dpn, []);    % dpn x Nnl
  Uds_qp = reshape(m.Qm*Und, m.dpn, []);  % dpn x Nnl

  [Fs_qp, Z, dFsdUs_qp, dFsdUds_qp, dFsdPars_qp] = m.fcont(Us_qp, Z, Uds_qp, Parsc);

  Fn = m.Tm*reshape(Fs_qp, [], 1);
  
  dFsdUs_qp = num2cell(dFsdUs_qp, [1, 2]);
  dFsdUds_qp = num2cell(dFsdUds_qp, [1, 2]);
  dFsdPars_qp = num2cell(dFsdPars_qp, [1, 2]);
  
  dFndUn   = m.Tm*blkdiag(dFsdUs_qp{:})*m.Qm;
  dFndUnd  = m.Tm*blkdiag(dFsdUds_qp{:})*m.Qm;
  dFndPars = m.Tm*blkdiag(dFsdPars_qp{:});

  if length(varargin)>=2  % Project/Mask everything back
    Fn = pU'*Fn;
    dFndUn = pU'*dFndUn*pU;
    dFndUnd = pU'*dFndUnd*pU;
    dFndPars = pU'*dFndPars;
  end
end
