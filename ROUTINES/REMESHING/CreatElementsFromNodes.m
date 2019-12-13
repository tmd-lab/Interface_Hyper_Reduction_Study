function [element_idx] = CreatElementsFromNodes(x,y)
% CreatElementsFromNodes: Creates Elements based on Node cloud
% INPUTS:
%   x : Vector containing all x - Coordinates
%   y : Vector containing all y - Coordinates
% OUTPUTS:
%   element_idx     : cell, containing Elements (Index of vertices)


[~,bolt_idx,pgon]=CreateBRBInterface(x,y);

%% Create elements on interface with Delaunay triangulation
triangle_idx=delaunay(x,y);   % Delaunay triangulation
% triplot(triangle_idx,BottomInterface.x,BottomInterface.y)

%Celan up triangles in the bolt holes
all_bolt_idx= cell2mat(bolt_idx);
delete_idx=logical(false(size(triangle_idx)));   % Initialize Matrix
for i=1:size(all_bolt_idx,1)
    delete_idx=delete_idx+(all_bolt_idx(i)==triangle_idx);
end
delete_triangle=all(delete_idx,2);
triangle_idx(delete_triangle==1,:)=[];


% Start algorithm to create patches
tic
element_idx{:,:}=[];
[pgon_red,element_idx,triangle] = RemoveTriangle(pgon,x,y,triangle_idx);
maxiter=25;  % Maximal number of iterations
i=1;
while ~isempty(triangle) && i<maxiter
    i=i+1;
    [pgon_red,element_idx,triangle] = RemoveTriangle(pgon_red,x,y,triangle,element_idx);
end
time=toc;
disp(['Time to Create a sqare mesh from Triangles ',num2str(time), 's']);
end

