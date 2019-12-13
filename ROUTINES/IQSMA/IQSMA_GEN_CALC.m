function [varargout] = IQSMA_GEN_CALC(pars, modelmats, expdat, Prestress, opts)
%IQSMA_CALC calculates the hysteresis loop and estimates the modal 
% characteristics.
%
% USAGE:
%   [Errors] = QSMA_CALC(pars, modelmats, expdat, Prestress, opts);
%       (or)
%   [Response, Errors] = QSMA_CALC(pars, modelmats, expdat, Prestress, opts);
% INPUTS:
%   pars            : (Npatches) x (Npars)
%   modelmats       : Structure
%                       M             : Mass matrix
%                       K             : Stiffness matrix
%                       L             : Null-space matrix
%                       Np            : Number of patches
%                       NLFORCINGFUNC : @(U, prev) Forcing function &
%                                       Jacobian
%                       CONTACTFUNC   : @(uxyn, uxyntxynp) Contact function
%                       NLRES         : @(NLFn, X, K, L, Fs, Fv, prev) 
%                                       Nonlinear residual function
%   expdat          : Structure
%   Prestress       : Scalar prestress value
%   opts            : Miscellaneous options for solution step
%                       No      : Number of quadrature points for hysteresis
%                                 integration
%                       mdid    : List of eigenmode numbers
%                       md      : Mode to choose from mdid list
%                       minlA   : Minimum Alpha (in logscale)
%                       maxlA   : Maximum Alpha (in logscale)
%                       Na      : Number of Alphas (in logscale)
%                       Display : 2 possible values tested
%                                'iter': display iteration information for
%                                           prestress and progress along
%                                           backbone
%                                'off' : no display
%                       LscPars  : List of parameter columns to be
%                                   interpreted as logscale (10^%T)
%                       SymL     : List for extended symmetry (see below
%                                   for example
%                       PConv    : Function handle for parameter basis
%                                   conversion. Empty if none necessary
%                                       
% OUTPUTS: 2 varargout modes:
%   nargout = 1; Use this for the errors (in log-scale) for Genetic algorithm
%       Errors = log10([Werr; Zerr])
%   nargout = 2; Use this for general response evaluation
%       Outputs = [Q W Z D];
%       Errors = [Werr; Zerr; Derr]
    %% Log-Scale Parameters
    pars(:, opts.LscPars) = 10.^(pars(:, opts.LscPars));
    if ~isempty(opts.PConv)
        pars = opts.PConv(pars);
    end
    %% Interface Parameters
    NLFUN = @(U, prev) modelmats.NLFORCINGFUNC(U, prev, pars);
    CFUNC = @(uxyn, uxyntxynp) modelmats.CONTACTFUNC(uxyn, uxyntxynp, pars);
    %% Prestress Analysis
    opt = optimoptions('fsolve', 'Display', opts.Display, 'SpecifyObjectiveGradient', true);
%     X0 = modelmats.K\(Prestress*modelmats.Fv);
    X0 = load('statsol.dat');
    prev.uxyntxyn = zeros(size(modelmats.QuadMats.Q,1)/3,6);
    
    [Xstat, ~, ~, ~, dR0] = fsolve(@(X) NLRES(NLFUN, [X; 0], modelmats.K, ...
        modelmats.L, Prestress*modelmats.Fv, modelmats.Fv*0, ...
        prev), X0, opt);
    [Px, Py, Pn] = CFUNC(reshape(modelmats.QuadMats.Q*L(1:(modelmats.MESH.dpn*modelmats.MESH.Nn),:)*Xstat, 3, [])', prev.uxyntxyn);
    Pstat = reshape([Px Py Pn]', size(modelmats.QuadMats.Q,1), 1);
    dlmwrite('statsol.dat', Xstat, 'precision', '%10e');
    
    prev.uxyntxyn = [reshape(modelmats.QuadMats.Q*L(1:(modelmats.MESH.dpn*modelmats.MESH.Nn),:)*Xstat, 3, [])' Px Py Pn];
    %% Zero-Amplitude Mode Shapes
    [Vm, Dm] = eigs(dR0, modelmats.M, 10, 'SM');
    [Dm, si] = sort(sqrt(abs(diag(Dm)))/(2*pi));
    V = Vm(:, opts.mdids(opts.md));
    V = V/(sqrt(V'*modelmats.M*V));
    %% QSMA
    Alphas = logspace(opts.minlA, opts.maxlA, opts.Na);
    
    Q = zeros(opts.Na, 1);
    W = Q;
    Z = Q;
    D = Q;
    Ro = zeros(size(modelmats.R,1), opts.Na);
    dx = zeros(size(modelmats.QuadMats.Q,1), opts.Na);
    dy = zeros(size(modelmats.QuadMats.Q,1), opts.Na);
    dn = zeros(size(modelmats.QuadMats.Q,1), opts.Na);
    
    X0 = Xstat + dR0\(Alphas(1)*modelmats.M*V);
    opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
    mshape = V;
    NLRESIDFUNC = @(U, alpha, prev) NLRES(NLFUN, [U; alpha], ...
            modelmats.K, modelmats.L, Prestress*modelmats.Fv, ...
            modelmats.M*mshape, prev);
%     parfor i=1:opts.Na
    for i=1:opts.Na
        [HYST, BB, ~, ~] = GENSTEPPEDHYSTERESIS(NLRESIDFUNC, CFUNC, ...
            Xstat+V*0, Alphas(i), opts.No, modelmats.MESH, modelmats.L, ...
            opt, modelmats.QuadMats,prev, 1);
        
        Ro(:, i) = modelmats.R*(BB.U(:,end)-Xstat);
        Q(i) = abs(mshape'*modelmats.M*(BB.U(:,end)-Xstat));
        W(i) = sqrt(2*Alphas(i)/range(mshape'*modelmats.M*(HYST.U-Xstat)));
        [D(i),dx(:,i),dy(:,i),dn(:,i)] = modelmats.DISSFUNC(modelmats.L(1:(modelmats.MESH.dpn*modelmats.MESH.Nn),:)*(HYST.U-Xstat), ...
            HYST.P-Pstat, modelmats.QuadMats, modelmats.MESH);
        Z(i) = D(i)/(2*pi*(Q(i)*W(i))^2);
        W(i) = W(i)/(2*pi);
        if ~strcmp(opts.Display, 'off')
            fprintf('Done %d/%d\n', i, opts.Na);
        end
    end
    
    %% Errors
    interped = interp1(Q, [W Z D], expdat.Q);
    Werr = rms(interped(:, 1)-expdat.W);
    Zerr = rms((interped(:, 2)-interped(1,2))-(expdat.Z-expdat.Z(1)));
    Derr = rms(interped(:, 3)-expdat.D);
    if sum(isnan([Werr Zerr Derr]))~=0
        disp('Isnan');
%         Werr = rms(expdat.W);
%         Zerr = rms(expdat.Z);
%         Derr = rms(expdat.D);
    end
    
    %% Outputs
    if nargout==1  %% Primary operation for Genetic Algorithm
        varargout{1} = log10([Werr;Zerr]);
    elseif nargout==2  %% Complete Simulation output
        varargout{1} = [Q W Z D];
        varargout{2} = [Werr; Zerr; Derr];
    elseif nargout==3  %% Complete Simulation output with largest hysteresis curve
        varargout{1} = [Q W Z D];
        varargout{2} = [Werr; Zerr; Derr];
        varargout{3} = Ro;
    elseif nargout==4  %% Complete Simulation output with largest hysteresis curve
        varargout{1} = [Q W Z D];
        varargout{2} = [Werr; Zerr; Derr];
        varargout{3} = Ro;        
        varargout{4} = struct('dx', dx, 'dy', dy, 'dn', dn);
    end
end