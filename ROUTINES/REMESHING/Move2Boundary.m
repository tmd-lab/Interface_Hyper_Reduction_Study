function [Moved_Nodes,Boundary_idx] = Move2Boundary(Nodes,Boundary,tolerance)
% Move2Boundary: Moves nodes close to boundary to boundary
% INPUTS:
%   Node:   Array with coordinates [x,y] 
%   Boundary : Corner nodes of Boundary [x,y]
%   tolerance: Tolerance for which nodes are moved
% OUTPUTS:
% New_Nodes: List of nodes with additionaly placed nodes


%   !! Important for a 4 Node input a boundary parallel to the coordinate
%   axis is assumed to reduce computational effort. 
%   Remove first if loop to make the code more generall.

Boundary=[Boundary;Boundary(1,:)];
Boundary_idx=[];

if length(Boundary)==4+1 
    for k=1:4
        idx=find(Boundary(k,:)==Boundary(k+1,:));
        if idx==1 || idx ==2
            idx_nodes=find(and(Boundary(k,idx)-tolerance<=Nodes(:,idx) , Boundary(k,idx)+ tolerance>=Nodes(:,idx) ));
            if ~isempty(idx_nodes)
                Nodes(idx_nodes,idx)=Boundary(k,idx);
                Boundary_idx=[Boundary_idx;idx_nodes];
            end
        else
            error('Your input is not consistent - Process terminated')
        end
    end
    Boundary_idx=unique(Boundary_idx);
else
    for k=1:length(Boundary)-1
        % Find potential points which are located close to the boundary
        dist=norm(Boundary(k,:)-Boundary(k+1,:));
        idx_candidates1= rangesearch(Nodes,Boundary(k,:),dist+tolerance);
        idx_candidates2= rangesearch(Nodes,Boundary(k+1,:),dist+tolerance);
        idx_candidate=intersect(cell2mat(idx_candidates1),cell2mat(idx_candidates2));
        for l=1:length(idx_candidate)
         % Calculating Distance to line between two boundary points
         a=Boundary(k,:); % First point of line
         n=(Boundary(k+1,:)-Boundary(k,:))/norm(Boundary(k+1,:)-Boundary(k,:)); % direction vector
         p=Nodes(idx_candidate(l),:); 
         project=(a-p)-(dot( (a-p), n ))*n; % Projection Vector
         newNode=p+project;
         
         if norm(Boundary(k,:)-Nodes(idx_candidate(l),:))<=tolerance
             Nodes(idx_candidate(l),:)=Boundary(k,:);
             Boundary_idx=[Boundary_idx;idx_candidate(l)];
             
         elseif norm(Boundary(k+1,:)-Nodes(idx_candidate(l),:))<=tolerance
             Nodes(idx_candidate(l),:)=Boundary(k+1,:);
             Boundary_idx=[Boundary_idx;idx_candidate(l)];
             
         elseif ( max(norm(Boundary(k,:)-newNode),norm(Boundary(k+1,:)-newNode))<=dist) ...
                 && norm(project)<= tolerance
             
             Nodes(idx_candidate(l),:)=newNode;
             Boundary_idx=[Boundary_idx;idx_candidate(l)];
         end
        end
        
    end
    Boundary_idx=unique(Boundary_idx);
    if length(Boundary_idx)>=3
    idx_sort=convhull(Nodes(Boundary_idx,1),Nodes(Boundary_idx,2));
    idx_sort=idx_sort(1:end-1);
    Boundary_idx=Boundary_idx(idx_sort);
    end
end
Moved_Nodes=Nodes;
end

