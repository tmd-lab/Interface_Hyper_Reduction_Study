function [Elements_idx] = FindInteriorNodes(Elements_idx,New_Nodes,Old_Nodes,toleranceCoincident)
% FindInteriorNodes: Findes nodes from lost that are inside a spesific
% element
% INPUTS:
%    Elements_idx: type cell (Number of Elements x 1): Each cell contins node ID of vertices
%    New_Nodes: (Number of New Nodes, 2): [x,y] Coord of new nodes 
%    Old_Nodes: (Number of Old Nodes, 2): [x,y] Coord of old nodes
%    toleranceCoincident: tolerance for points are considered coincident
% OUTPUTS:
%    Elements_idx: typ cell (Number of Elements x 2): the second
%    column of each cell contains a vector with ID of old nodes located inside
%    new elements

for  ii=1:length(Elements_idx)

    
    VerticesNodes=[New_Nodes(Elements_idx{ii,1},1),New_Nodes(Elements_idx{ii,1},2)];
    if length(Elements_idx{ii,1})==4 % Quad
        
        %Find nodes which are coincident
        DiffC1=Old_Nodes-VerticesNodes(1,:);
        DiffC2=Old_Nodes-VerticesNodes(2,:);
        DiffC3=Old_Nodes-VerticesNodes(3,:);
        DiffC4=Old_Nodes-VerticesNodes(4,:);
        
        ison=(vecnorm(DiffC1,2,2)<toleranceCoincident | vecnorm(DiffC2,2,2)<toleranceCoincident | ...
            vecnorm(DiffC3,2,2)<toleranceCoincident | vecnorm(DiffC4,2,2)<toleranceCoincident);
        
        direction1=VerticesNodes(end,:)-VerticesNodes(1,:);
        direction2=VerticesNodes(2,:)-VerticesNodes(1,:);
        direction3=VerticesNodes(2,:)-VerticesNodes(3,:);
        direction4=VerticesNodes(end,:)-VerticesNodes(3,:);
        
        angle1=atan2(direction1(2),direction1(1));
        angle2=atan2(direction2(2),direction2(1));
        angle3=atan2(direction3(2),direction3(1));
        angle4=atan2(direction4(2),direction4(1));
        
        
        nodesangle1=atan2(DiffC1(:,2),DiffC1(:,1));
        nodesangle3=atan2(DiffC3(:,2),DiffC3(:,1));
        
        if abs(angle1-angle2)<pi
            criteria1=(nodesangle1<=angle1 & nodesangle1>=angle2);
        else
            if angle1<0
                angle1=angle1+2*pi;
            end
            if angle2<0
                angle2=angle2+2*pi;
            end
            nodesangle1(nodesangle1<0)=nodesangle1(nodesangle1<0)+2*pi;
            
            criteria1=(nodesangle1<=angle1 & nodesangle1>=angle2);
            
        end
        
        if abs(angle3-angle4)<pi
            criteria2=(nodesangle3<=angle3 & nodesangle3>=angle4);
        else
            if angle3<0
                angle3=angle3+2*pi;
            end
            if angle4<0
                angle4=angle4+2*pi;
            end
            nodesangle3(nodesangle3<0)=nodesangle3(nodesangle3<0)+2*pi;
            
            criteria2=(nodesangle3<=angle3 & nodesangle3>=angle4);
        end
        Elements_idx{ii,2}=find(or(ison, and(criteria1 , criteria2)));
        

    elseif length(Elements_idx{ii,1})==3   
        %Find nodes which are coincident
        DiffC1=Old_Nodes-VerticesNodes(1,:);
        DiffC2=Old_Nodes-VerticesNodes(2,:);
        DiffC3=Old_Nodes-VerticesNodes(3,:);
        
        ison=(vecnorm(DiffC1,2,2)<toleranceCoincident | vecnorm(DiffC2,2,2)<toleranceCoincident | ...
            vecnorm(DiffC3,2,2)<toleranceCoincident);
        
        direction1=VerticesNodes(end,:)-VerticesNodes(1,:);
        direction2=VerticesNodes(2,:)-VerticesNodes(1,:);
        direction3=VerticesNodes(2,:)-VerticesNodes(end,:);
        direction4=VerticesNodes(1,:)-VerticesNodes(end,:);
        
        angle1=atan2(direction1(2),direction1(1));
        angle2=atan2(direction2(2),direction2(1));
        angle3=atan2(direction3(2),direction3(1));
        angle4=atan2(direction4(2),direction4(1));
        
        
        nodesangle1=atan2(DiffC1(:,2),DiffC1(:,1));
        nodesangle3=atan2(DiffC3(:,2),DiffC3(:,1));
        
        if abs(angle1-angle2)<pi
            criteria1=(nodesangle1<=angle1 & nodesangle1>=angle2);
        else
            if angle1<0
                angle1=angle1+2*pi;
            end
            if angle2<0
                angle2=angle2+2*pi;
            end
            nodesangle1(nodesangle1<0)=nodesangle1(nodesangle1<0)+2*pi;
            
            criteria1=(nodesangle1<=angle1 & nodesangle1>=angle2);
            
        end
        
        if abs(angle3-angle4)<pi
            criteria2=(nodesangle3<=angle3 & nodesangle3>=angle4);
        else
            if angle3<0
                angle3=angle3+2*pi;
            end
            if angle4<0
                angle4=angle4+2*pi;
            end
            nodesangle3(nodesangle3<0)=nodesangle3(nodesangle3<0)+2*pi;
            
            criteria2=(nodesangle3<=angle3 & nodesangle3>=angle4);
        end
        Elements_idx{ii,2}=find(or(ison, and(criteria1 , criteria2)));
    else
        error('');
    end    

end

end

