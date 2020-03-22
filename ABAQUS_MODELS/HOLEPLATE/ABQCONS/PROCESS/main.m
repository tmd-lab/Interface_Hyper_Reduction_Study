clc
clear all
addpath('../../../../ROUTINES/')
addpath('../../../../ROUTINES/FEM/')

Mdat = dlmread('../MATEX/MATEX_MASS1.mtx');
Kdat = dlmread('../MATEX/MATEX_STIF1.mtx');
Nn = max([max(max(Mdat(:, [1 3]))), max(max(Kdat(:, [1 3])))]);
Nlag = abs(min([min(min(Mdat(:, [1 3]))), min(min(Kdat(:, [1 3])))]));

Ndi = dlmread('../MATEX/NodeLabels.dat');

Nd_dfs = zeros(Nn+Nlag, 2);  % Node labels and corresponding number of DoFs
Nd_dfs(:, 1) = unique([reshape(Mdat(:, [1 3]), [],1); reshape(Kdat(:, [1 3]), [],1)]);
for i=1:length(Nd_dfs)
    Nd_dfs(i,2) = max(Mdat(Mdat(:, 1)==Nd_dfs(i,1), 2));
end
Nd_dfs = [Nd_dfs(Nlag+1:end, :); Nd_dfs(Nlag:-1:1, :)];

%% Construct Matrices
K = zeros((Nn-Nlag)*3+2*Nlag*6, (Nn-Nlag)*3+2*Nlag*6, 'single');
M = zeros((Nn-Nlag)*3+2*Nlag*6, (Nn-Nlag)*3+2*Nlag*6, 'single');
for n=1:size(Kdat, 1)
    ni = find(Kdat(n,1)==Nd_dfs(:,1));
    nj = find(Kdat(n,3)==Nd_dfs(:,1));
    
    i = sum(Nd_dfs(1:(ni-1), 2))+Kdat(n,2);
    j = sum(Nd_dfs(1:(nj-1), 2))+Kdat(n,4);
    
    K(i, j) = Kdat(n, 5);
    K(j, i) = Kdat(n, 5);
end
for n=1:size(Mdat, 1)
    ni = find(Mdat(n,1)==Nd_dfs(:,1));
    nj = find(Mdat(n,3)==Nd_dfs(:,1));
    
    i = sum(Nd_dfs(1:(ni-1), 2))+Mdat(n,2);
    j = sum(Nd_dfs(1:(nj-1), 2))+Mdat(n,4);
    
    M(i, j) = Mdat(n, 5);
    M(j, i) = Mdat(n, 5);
end

disp('Done!')
K = sparse(double(K));
M = sparse(double(M));
%% Reorder
T = eye(Nn-Nlag);
Ndo = setdiff((1:Nn-Nlag)', Ndi);
% T = T([Ndi; Ndo], :);
T = T(:, [Ndi; Ndo]);
T = blkdiag(kron(T, eye(3)), eye(Nlag*2*6));
T = sparse(T);

K = T'*K*T;
M = T'*M*T;

%% Plot structure
or = load('../../PROCESS/EXTRACTION.mat', 'K', 'M');

figure(1)
spy(blkdiag(or.K, zeros(2*Nlag*6)))
figure(2)
spy(K)
figure(3)
spy(K-blkdiag(or.K, zeros(2*Nlag*6)))
figure(4)
spy(K([(1:length(Ndi)*3) end-2*Nlag*6+1:end],[(1:length(Ndi)*3) end-2*Nlag*6+1:end]))
figure(5)
spy(K([(1:length(Ndi)*3) end-2*Nlag*6+1:end],[(1:length(Ndi)*3) end-2*Nlag*6+1:end])-blkdiag(or.K((1:length(Ndi)*3),(1:length(Ndi)*3)), zeros(2*Nlag*6)))
