clc
clear all
addpath('../../../ROUTINES/')
addpath('../../../ROUTINES/FEM/')

Mdat = dlmread('../MATEX/MATEX_MASS1.mtx');
Kdat = dlmread('../MATEX/MATEX_STIF1.mtx');
Nn = max([max(max(Mdat(:, 1:2))), max(max(Kdat(:, 1:2)))]);
Ndi = dlmread('../MATEX/NodeLabels.dat');

%% Construct Matrices
K = zeros(Nn*3, Nn*3, 'single');
M = zeros(Nn*3, Nn*3, 'single');
for n=1:size(Kdat,1)
    i = (Kdat(n,1)-1)*3+Kdat(n,2);
    j = (Kdat(n,3)-1)*3+Kdat(n,4);
    
    K(i, j) = Kdat(n, 5);
    K(j, i) = Kdat(n, 5);
end
for n=1:size(Mdat,1)
    i = (Mdat(n,1)-1)*3+Mdat(n,2);
    j = (Mdat(n,3)-1)*3+Mdat(n,4);
    5
	M(i, j) = Mdat(n, 5);
    M(j, i) = Mdat(n, 5);
end
disp('Done!')
K = sparse(double(K));
M = sparse(double(M));

%% Reorder
T = eye(Nn);
Ndo = setdiff((1:Nn)', Ndi);
% T = T([Ndi; Ndo], :);
T = T(:, [Ndi; Ndo]);
T = sparse(kron(T, eye(3)));

K = T'*K*T;
M = T'*M*T;

%% Construct Mesh Structure
Els = dlmread('../MATEX/Elements.dat');
MESH.Nds  = dlmread('../MATEX/Nodes.dat');
MESH.Quad = Els(:,1:end-1)+1;
MESH.Tri  = [];
MESH.Nn      = size(MESH.Nds, 1);
MESH.Ne_Quad = size(MESH.Quad, 1);
MESH.Ne_Tri  = size(MESH.Tri, 1);
MESH.Ne      = MESH.Ne_Quad + MESH.Ne_Tri;

Pels = cell(5, 1);
Pnds = cell(5, 1);
figure(1);
clf()
for i = 1:5
    Pels{i} = find(Els(:, end)==i-1);
    Pnds{i} = unique(MESH.Quad(Pels{i}, 2:end));
    
    SHOW2DMESH(MESH.Nds, MESH.Tri, [(1:length(Pels{i}))' MESH.Quad(Pels{i},2:end)], i, -1, -100, MESH.Ne)
end
axis equal
axis off

save('EXTRACTION.mat', 'MESH', 'Pels', 'Pnds', 'K', 'M', 'Ndi', 'T')

%% Check Mode shapes
[V, D] = eigs(K, M, 30, 'SM');
V = V./sqrt(diag(V'*M*V)');
%%
mi = 7; figure(2); clf(); SHOW3D(MESH.Nds, MESH.Tri, MESH.Quad, V(1:3:MESH.Nn*3,mi), V(2:3:MESH.Nn*3,mi), V(3:3:MESH.Nn*3,mi), 0.001, V(3:3:MESH.Nn*3,mi)); axis equal; title(sprintf('Mode %d: W = %f Hz', mi, sqrt(D(mi,mi))/2/pi));