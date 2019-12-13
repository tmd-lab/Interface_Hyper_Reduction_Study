function [ReducedNodes] = RemoveCloseNodes(Nodes,tolerance)
% RemoveCloseNodes:If two nodes are to close togher on will be removed
% INPUTS:
%   Nodes: Array with nodes [x,y]
%   tolerance: Tolerance when nodes should be combined
% OUTPUTS:
%   ReducedNodes : Reduced array of nodes


Dist=ones(length(Nodes))*Inf;

for k=1:length(Nodes)
    for l=k+1:length(Nodes)
        Dist (l,k)=norm(Nodes(k,:)-Nodes(l,:));
    end
end
[idx_col,idx_row]=find(Dist<=tolerance);
if length([idx_col;idx_row])~=length(unique([idx_col;idx_row]))
    error('Needs to be implemented in case this case ocurrs.');
else
    idx_remove=[idx_col,idx_row];
    for k=1:size(idx_remove,1)
        x_new=mean(Nodes(idx_remove(k,:),1));
        y_new=mean(Nodes(idx_remove(k,:),2));
        Nodes(idx_remove(k,1),:)=[x_new ,y_new];
        Nodes(idx_remove(k,2),:)=[Inf, 0];
    end
    idx_del=Nodes(:,1)==Inf;
    Nodes(idx_del,:)=[];
end

ReducedNodes=Nodes;
end


