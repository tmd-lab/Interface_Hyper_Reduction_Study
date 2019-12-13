function [R,dRdU,dRda,dRdX] = NLRES_RED(nlforcingfunc, U, K, L, Lt, Fl, prev)
%NLRES_RED Returns the nonlinear residue and its Jacobians for a
%quasi-static simulation. The parameter by default is the
%scaling applied to the given forcing vector.
% USAGE:
%	[R,dRdU,dRda] = NLRES(contactfunc, U, K, L, Fl, prev); 
% INPUTS:
%   contactfunc	:
%   U		: ((Ndof+1)x1) displacement vector [U; alpha]
%   K		: (NdofxNdof) stiffness matrix
%   L		: (NphxNdof) Null transformation matrix
%   Fl		: (Ndofx1) force vector (to be scaled by alpha)
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

    Uph = L*(U(1:end-1));
    [Fnl, Jnl] = nlforcingfunc(Uph,prev);

    R 		= K*(U(1:end-1)) + Lt'*Fnl - U(end)*Fl;
    dRdU 	= K + Lt'*Jnl*L;
    dRda 	= -Fl;
    dRdX 	= [dRdU dRda];
end