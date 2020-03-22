clc
clear all
addpath('../../../ROUTINES/')
addpath('../../../ROUTINES/FEM/')

load('EXTRACTION.mat', 'MESH', 'Pels', 'Pnds', 'K', 'M', 'Ndi', 'T')
Nn = size(K,1)/3;
Np = length(Pels);
[Q1, T1] = ZTE_ND2QP(MESH, 1);

%% Patch Reduction
Area = cell(1, Np);
qps = cell(1, Np);
ctrds = zeros(Np, 3);

colormap(jet(Np));
colos = colormap;

figure(1)
clf()
for n=1:Np
    Area{n} = T1(Pnds{n}, Pels{n})*ones(length(Pels{n}), 1);
    qps{n} = Q1(Pels{n}, :)*MESH.Nds;
    ctrds(n, :) = sum(Area{n}.*MESH.Nds(Pnds{n}, :))/sum(Area{n});
    
    SHOW2DMESH(MESH.Nds, zeros(0,4), MESH.Quad(Pels{n},:), n, -1, -100)
    plot(ctrds(n, 1), ctrds(n, 2), 'ko', 'MarkerFaceColor', colos(n,:))
end
PatchAreas = cellfun(@(c) sum(c), Area);
axis equal; axis off; colormap(jet(Np))

%% Constraint Weak Form Matrices
NTN = sparse(MESH.Nn*3, MESH.Nn*3);
NTG = sparse(MESH.Nn*3, Np*6);
GTG = sparse(Np*6, Np*6);
BNV = sparse(MESH.Nn, Np);  % Force Transformation Matrix
for n=1:Np
    [P, Nums, NTNmat, NTGmat, GTGmat] = CONSPATCHMAT(MESH.Nds, [], MESH.Quad(Pels{n}, :), ctrds(n, :));
    
    NTN = NTN + NTNmat;
    NTG(:, (n-1)*6+(1:6)) = NTG(:, (n-1)*6+(1:6)) + NTGmat;
    GTG((n-1)*6+(1:6), (n-1)*6+(1:6)) = GTG((n-1)*6+(1:6), (n-1)*6+(1:6)) + GTGmat;
    
    BNV(Pnds{n}, n) = 1.0/sum(Area{n});
end
NTN = 0.5*(NTN+NTN');
%% ROM Development
cnum = 1e12;

M1 = sparse(blkdiag(M, zeros(2*Np*6)));
K1 = sparse([K(1:MESH.Nn*3,:) zeros(MESH.Nn*3, Np*6) -cnum*NTN;
      K(MESH.Nn*3+1:end,:) zeros((Nn-MESH.Nn)*3,Np*6+MESH.Nn*3);
      zeros(Np*6, Np*6+Nn*3) cnum*NTG';
      -cnum*NTN' zeros(MESH.Nn*3, (Nn-MESH.Nn)*3) cnum*NTG zeros(MESH.Nn*3)]);
disp(cond(K1))
% HCB Here
Ngen = 20;
[Mh, Kh, Th] = HCBREDUCE(M1, K1, Nn*3+(1:Np*6), Ngen);
disp('HCB Done')
%% Compare Modes
[V, D] = eigs(K, M, size(Kh,1), 'SM');
% [V, D] = eig(full(K), full(M));
[D, si] = sort(diag(D));
V = V(:, si);
V = V./sqrt(diag(V'*M*V)');

[Vh, Dh] = eig(full(Kh), full(Mh));
[Dh, si] = sort(diag(Dh));
Vh = Vh(:, si);
Vh = Vh./sqrt(diag(Vh'*Mh*Vh)');

mhs = find(isfinite(Dh));
VH = Th(1:Nn*3,:)*Vh;

Nv = length(mhs)+5;
figure(2)
br=bar3((VH(:, mhs)'*M*V(:, 1:Nv)).^2./(diag(VH(:,mhs)'*M*VH(:,mhs)).*diag(V(:, 1:Nv)'*M*V(:, 1:Nv))'));
for k=1:size(br,2); br(k).CData = br(k).ZData; br(k).FaceColor = 'interp'; end
xlabel('Original Model')
ylabel('ROM')
xlim([1 Nv])
ylim([1 Nv])

%% Plot modes

mi = 9; 
figure(3); clf(); SHOW3D(MESH.Nds, MESH.Tri, MESH.Quad, V(1:3:MESH.Nn*3,mi), V(2:3:MESH.Nn*3,mi), V(3:3:MESH.Nn*3,mi), 0.001, V(3:3:MESH.Nn*3,mi)); axis equal; title(sprintf('Mode %d: W = %f Hz', mi, sqrt(D(mi))/2/pi)); colormap(jet(11))
figure(4); clf(); SHOW3D(MESH.Nds, MESH.Tri, MESH.Quad, real(VH(1:3:MESH.Nn*3,mi)), real(VH(2:3:MESH.Nn*3,mi)), real(VH(3:3:MESH.Nn*3,mi)), 0.001, real(VH(3:3:MESH.Nn*3,mi))); axis equal; title(sprintf('Mode %d: W = %f Hz', mi, sqrt(Dh(mi))/2/pi)); colormap(jet(11))