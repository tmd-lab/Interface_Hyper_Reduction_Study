clc
clear all
addpath('../ROUTINES/ROUTINES')
addpath('../ROUTINES/ROUTINES/FEM')
addpath('../ROUTINES/ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/ROUTINES/QUASISTATIC')
addpath('../ROUTINES/ROUTINES/TRANSIENT')
addpath('../ROUTINES/ROUTINES/SOLVERS')

% ROM levels
sel_method = 'P';
Nels = 52;  % [52 100 140 204 240 304]

% sel_method = 'U';
% Nels = 568; % [36 74 122 232 448 568 588]
% for Nels=[36 74 122 232 448 568 588]

% sel_method = 'PD';
% Nels = 252;  % [48 108 152 200 252 292]

%% Load Matrices and setup MESH
load(sprintf('../FULLJOINT_ROM_PREPARE/ROMS/ROM_%s_%dELS.mat',sel_method,Nels), 'M', 'K', 'R', 'L', 'Fv', 'TFMf', 'TFMfhcb', 'MESH');
MESH = MESH2D(MESH.Nds, 3, [], MESH.Quad, 2);
MESH = MESH.SETCFUN(@(u, z, ud, P) ELDRYFRICT(u, z, ud, P, 0), sparse(2, MESH.Ne*MESH.Nq^2));  % Contact Function
MESH.Ltran = TFMfhcb;