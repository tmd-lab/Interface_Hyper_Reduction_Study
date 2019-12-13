function [R,dRdUl,dRdq,dRdUlq] = QSHMARES(nlforcingfunc, Ulq, M, K, L, Fs, Us, prev, varargin)
%QSHMARES Returns the nonlinear residue and its Jacobians for a
%quasi-static hysteretic modal analysis simulation. The parameter by 
%default is the
%scaling applied to the given forcing vector.
% USAGE:
%	[R,dRdUw,dRdq] = QSHMARES(contactfunc, Uwq, K, L, Fl, prev, Cmat); 
% INPUTS:
%   contactfunc	:
%   Uwq		: ((Ndof+2)x1) displacement vector with frequency and amplitude [U; w; q]
%   K		: (NdofxNdof) stiffness matrix
%   L		: (NphxNdof) Null transformation matrix
%   Fs		: (Ndofx1) force vector (to be held constant)
%   Us		: (Ndofx1) static solution vector
%   prev        : Previous state structure
%   	ux      : (Nnx1) vectors of interface nodal rel. x displacements
%	uy      : (Nnx1) vectors of interface nodal rel. y displacements  
%	un      : (Nnx1) vectors of interface nodal rel. normal displacements
%   	tx	: (Nnx1) vectors of interface nodal rel. x tractions
%       ty	: (Nnx1) vectors of interface nodal rel. y tractions
%       tn      : (Nnx1) vectors of interface nodal rel. normal tractions
% OUTPUTS:
%   R		: (Ndofx1) vector of residue
%   dRdX 	: (Ndofx(Ndof+1)) matrix of combined Jacobians    
%   dRdU	: (NdofxNdof) matrix of derivatives of R w.r.t
%   			displacements
%   dRda	: (Ndofx1) vector of derivative of R w.r.t alpha
    if nargin>=9
        Ulq = varargin{1}.*Ulq;
    end
    mspc = 1.0;
    if nargin>=10
        mspc = varargin{2};
    end
    U = Ulq(1:end-2); lam = Ulq(end-1);
    if nargin>=11
        q = 10^Ulq(end);
        dqdu = q*log(10);
    else
        q = Ulq(end);
        dqdu = 1.0;
    end

    Uph = L*U;
    [Fnl, Jnl] = nlforcingfunc(Uph,prev);

    R 		= [K*U+L'*Fnl-lam*M*(U-Us)-Fs; 
               mspc*((U-Us)'*M*(U-Us)-q^2)];
    dRdUl 	= [K+L'*Jnl*L-lam*M, -M*(U-Us);
               2*mspc*(U-Us)'*M, 0];
    dRdq 	= [zeros(size(U)); -2*mspc*dqdu*q];
    if nargin>=9
        dRdUl = dRdUl.*varargin{1}(1:end-1);
        dRdq  = dRdq*varargin{1}(end);
    end

    dRdUlq 	= [dRdUl dRdq];
end