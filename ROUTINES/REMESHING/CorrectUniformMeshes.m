function [] = CorrectUniformMeshes(Nels)
% CorrectUniformMeshes: Correction of Node vector for unoform mesh.
% opens and deletes corresponding file
% INPUTS:
%   Nels : No of Elements

Ntyps = length(Nels);

for i=1:Ntyps
    load(sprintf('UniformMesh%dElements.mat',Nels(i)));
    
    Nds = OutputNodes;
    Quad_Els = OutputQuad;
    Tri_Els = OutputTri;
    A = FEMGRAPHADJ(Nds(:,2:3),Quad_Els,Tri_Els);
    G = graph(A);
    [bins,binsizes] = conncomp(G);
    reqdnds = find(bins==find(binsizes~=1));
    nreqdnds = setdiff(1:length(Nds),reqdnds);
    
    naa = zeros(size(Nds,1),1);
    naa(reqdnds) = 1:length(reqdnds);
    Nds(nreqdnds,:) = [];
    Quad_Els(:,2:5) = naa(Quad_Els(:,2:5));
    Tri_Els(:,2:4) = naa(Tri_Els(:,2:4));
    
    OutputNodes = Nds;
    OutputQuad = Quad_Els;
    OutputTri = Tri_Els;
    
    save(sprintf('UniformMesh%dElements.mat',Nels(i)),'OutputNodes','OutputQuad','OutputTri','NodeCorrelation');

end

