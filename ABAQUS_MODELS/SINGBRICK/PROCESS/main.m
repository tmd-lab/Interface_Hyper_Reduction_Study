clc
clear all

%% Construct Matrices - Unconstrained Case
Mdat = dlmread('../Job-1_MASS2.mtx');
Kdat = dlmread('../Job-1_STIF2.mtx');

Nn = max([max(max(Mdat(:, [1 3]))), max(max(Kdat(:, [1 3])))]);
Nlag = 0;

Nd_dfs = zeros(Nn+Nlag, 2);  % Node labels and corresponding number of DoFs
Nd_dfs(:, 1) = unique([reshape(Mdat(:, [1 3]), [],1); reshape(Kdat(:, [1 3]), [],1)]);
for i=1:length(Nd_dfs)
    Nd_dfs(i,2) = max(Mdat(Mdat(:, 1)==Nd_dfs(i,1), 2));
end
Nd_dfs = [Nd_dfs(Nlag+1:end, :); Nd_dfs(Nlag:-1:1, :)];

K = zeros(Nn*3, Nn*3, 'single');
M = zeros(Nn*3, Nn*3, 'single');
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

% Reorder
Ndis = [1 2 5 6];
T = eye(8);
T = T(:, [Ndis setdiff(1:8, Ndis)]);
T = kron(T, eye(3));

K = T'*K*T;
M = T'*M*T;

%% Construct Matrices - Constrained Case
Mdat = dlmread('../Job-2_MASS2.mtx');
Kdat = dlmread('../Job-2_STIF2.mtx');

Nn = max([max(max(Mdat(:, [1 3]))), max(max(Kdat(:, [1 3])))]);
Nlag = abs(min([min(min(Mdat(:, [1 3]))), min(min(Kdat(:, [1 3])))]));

Nd_dfs = zeros(Nn+Nlag, 2);  % Node labels and corresponding number of DoFs
Nd_dfs(:, 1) = unique([reshape(Mdat(:, [1 3]), [],1); reshape(Kdat(:, [1 3]), [],1)]);
for i=1:length(Nd_dfs)
    Nd_dfs(i,2) = max(Mdat(Mdat(:, 1)==Nd_dfs(i,1), 2));
end
Nd_dfs = [Nd_dfs(Nlag+1:end, :); Nd_dfs(Nlag:-1:1, :)];

Kc = zeros(Nn*3, Nn*3, 'single');
Mc = zeros(Nn*3, Nn*3, 'single');
for n=1:size(Kdat, 1)
    ni = find(Kdat(n,1)==Nd_dfs(:,1));
    nj = find(Kdat(n,3)==Nd_dfs(:,1));
    
    i = sum(Nd_dfs(1:(ni-1), 2))+Kdat(n,2);
    j = sum(Nd_dfs(1:(nj-1), 2))+Kdat(n,4);
    
    Kc(i, j) = Kdat(n, 5);
    Kc(j, i) = Kdat(n, 5);
end
for n=1:size(Mdat, 1)
    ni = find(Mdat(n,1)==Nd_dfs(:,1));
    nj = find(Mdat(n,3)==Nd_dfs(:,1));
    
    i = sum(Nd_dfs(1:(ni-1), 2))+Mdat(n,2);
    j = sum(Nd_dfs(1:(nj-1), 2))+Mdat(n,4);
    
    Mc(i, j) = Mdat(n, 5);
    Mc(j, i) = Mdat(n, 5);
end

disp('Done!')
Kc = sparse(double(Kc));
Mc = sparse(double(Mc));

Kc(1:8*3, :) = T'*Kc(1:8*3,:);
Mc(1:8*3, :) = T'*Mc(1:8*3,:);
Kc(:, 1:8*3) = Kc(:, 1:8*3)*T;
Mc(:, 1:8*3) = Mc(:, 1:8*3)*T;

save('MATS.mat', 'Kc', 'Mc', 'K', 'M', 'T', 'Ndis')
%% Plot structure
Ktmp = blkdiag(K, zeros(2*Nlag*6));
figure(10)
spy(Ktmp)
figure(20)
spy(Kc)
figure(30)
spy(Kc-Ktmp)
