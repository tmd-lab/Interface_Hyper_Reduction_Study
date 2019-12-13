function [MirroredNodes,MirroredElements] = MirrorElements(Nodes,Elements,BoundNodesidx,Direction,Offset)
% MirrorElements: Generates full model from quarter model
% INPUTS:
%   Nodes : Array containing nodes
%   Elements: Cell containing all elements	
%   BoundNodesidx: Index of nodes that are located on mirror axis
%   Direction: Mirror Axis ['X','Y']
%   Offset: Offset of mirror axis to coordinate axis
% OUTPUTS:
%   MirroredNodes : Array containing all nodes 
%   MirroredElements: cell containing all elements

N=size(Nodes,1);
NEle=length(Elements);
%Mirror  Node Vector
if Direction =='Y'
    New_Nodes=[Nodes;[Nodes(:,1),(-(Nodes(:,2)-Offset))+Offset]]; % y-Direction
elseif Direction=='X'
    New_Nodes=[Nodes;[(-(Nodes(:,1)-Offset))+Offset,Nodes(:,2)]]; % x-Direction
else
    error('Invalid input');
end

% Mirror Elements
for k=1:NEle
    idx_int=intersect(Elements{k,1},BoundNodesidx);
    if isempty(idx_int)
        Elements{k+NEle,1}= Elements{k,1}+N;
    else
        Elements{k+NEle,1}= Elements{k,1};
        Elements{k+NEle,1}(all(Elements{k,1}~=idx_int,1))= Elements{k,1}(all(Elements{k,1}~=idx_int,1))+N;
    end
end

%% Postprocessing
% Remove unused nodes from Node Matrix 
New_Nodes(BoundNodesidx+N,:)=[];

% Sort boundary (important for ID update)
BoundNodesidx=sort(BoundNodesidx,'descend'); 

% Update Element IDs
for k=NEle+1:length(Elements)
    for l=1:length(BoundNodesidx)
        idx=find(Elements{k,1} >= BoundNodesidx(l)+N);
      Elements{k,1} (idx )= Elements{k,1} (idx )-1;
    end
end

% Sort Elements Counter clockwise 
for k=1:length(Elements)
    if ispolycw(New_Nodes(Elements{k,1},1)',New_Nodes(Elements{k,1},2)')  
        Elements{k,1}= flip (Elements{k,1});
    end
end


MirroredNodes=New_Nodes;
MirroredElements=Elements;
end

