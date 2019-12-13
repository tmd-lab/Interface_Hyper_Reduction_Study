function [path,element_path] = Node2NodePath(Method,start_node,end_node,element,x,y)
% Node2NodePath:
% This funcion returns the vector with all nodes on the way between two
% nodes for a mesh (Method: dist)
% This function returns the vector with all nodes on the way from one node
% to the boundary  ( Method: path)
% INPUTS:
%   start_node: Index of start node
%   end_node: index of end node	
%   element: cell array containing all elements
%   x: Vector of x Coordinates 
%   y: Vector of y Coordinates 
% OUTPUTS:
%   path : node list / array containing all nodes along the path
%   element_path: list with elements on the path (only unique for along a
%   boundary)


element_path=[];
path=[start_node ];
if ~isempty(intersect(start_node,end_node))
   error('Invalid Input')
end

% Find Element of fist node
for i=1:size(element,1)
    idx(i)=length(intersect(start_node,element{i,1}))>=1;
end
idx=find(idx);


if strcmp(Method,'path')
    dx=diff(x([start_node,end_node]));
    dy=diff(y([start_node,end_node]));
    angle1=atan(dy/dx);
    angle2=atan2(dy , dx);

    distance=zeros(length(idx)*2,4);
    for p=1:length(idx)
        iidx(p)=find(element{idx(p),1}==start_node);
        if iidx(p)==length(element{idx(p),1})
            distance(2*p-1,1)=x(element{idx(p),1}(iidx(p)-1))-x(start_node);
            distance(2*p-1,2)=y(element{idx(p),1}(iidx(p)-1))-y(start_node);
            distance(2*p,1)=x(element{idx(p),1}(1))-x(start_node);
            distance(2*p,2)=y(element{idx(p),1}(1))-y(start_node);
        elseif iidx(p)==1
            distance(2*p-1,1)=x(element{idx(p),1}(end))-x(start_node);
            distance(2*p-1,2)=y(element{idx(p),1}(end))-y(start_node);
            distance(2*p,1)=x(element{idx(p),1}(iidx(p)+1))-x(start_node);
            distance(2*p,2)=y(element{idx(p),1}(iidx(p)+1))-y(start_node);
        else
            distance(2*p-1,1)=x(element{idx(p),1}(iidx(p)-1))-x(start_node);
            distance(2*p-1,2)=y(element{idx(p),1}(iidx(p)-1))-y(start_node);
            distance(2*p,1)=x(element{idx(p),1}(iidx(p)+1))-x(start_node);
            distance(2*p,2)=y(element{idx(p),1}(iidx(p)+1))-y(start_node);
        end
        b1=[distance(2*p-1,1) distance(2*p-1,2)];
        b2=[distance(2*p,1) distance(2*p,2)];
        % Check for same orientation
        distance(2*p-1,3)=abs(atan(b1(2)/b1(1))-angle1);
        if distance(2*p-1,3)>pi()/2
            distance(2*p-1,3)= pi()-distance(2*p-1,3);
        end
        distance(2*p,3)=abs(atan(b2(2)/b2(1))-angle1);
        if distance(2*p,3)>pi()/2
            distance(2*p,3)= pi()-distance(2*p,3);
        end
        %check for same direction
        distance(2*p-1,4)=atan2(b1(2),b1(1))-angle2;
        distance(2*p,4)=atan2(b2(2),b2(1))-angle2;

    end
    idx_temp=find(distance(:,3)<=pi()/4);
    [~,idx_temp2]=min(abs(distance(idx_temp,4)));
    idx_min=idx_temp(idx_temp2);
    r=round(idx_min/2);
    
elseif strcmp(Method,'all')
    r=1:length(idx);
    for j=1:length(idx)
        iidx(j)=find(element{idx(j),1}==start_node);
    end
    idx_min=1;

else
    error('Selected Method does not exist')
end

for k=1:length(r)
    p=r(k);
    element_path=[idx(p)];

    if mod(idx_min,2)==1
        direction=1;
        if iidx(p)== 1
            connection_nodes=[element{idx(p),1}(end) element{idx(p),1}(end-1)];
            path=[path element{idx(p),1}(end)];
        elseif iidx(p)== 2
            connection_nodes=[element{idx(p),1}(1) element{idx(p),1}(end)];
            path=[path element{idx(p),1}(iidx(p)-1)];
        else
            connection_nodes=[element{idx(p),1}(iidx(p)-1) element{idx(p),1}(iidx(p)-2)];
            path=[path element{idx(p),1}(iidx(p)-1)];
        end
        
    elseif mod(idx_min,2)==0
        direction=2;
        if iidx(p)== length(element{idx(p),1})
            connection_nodes=[element{idx(p),1}(1) element{idx(p),1}(2)];
            path=[path element{idx(p),1}(1)];
        elseif iidx(p)== length(element{idx(p),1})-1
            connection_nodes=[element{idx(p),1}(end) element{idx(p),1}(1)];
            path=[path element{idx(p),1}(end)];
        else
            connection_nodes=[element{idx(p),1}(iidx(p)+1) element{idx(p),1}(iidx(p)+2)];
            path=[path element{idx(p),1}(iidx(p)+1)];
        end
    else
        warning('Something is wrong');
    end
    
    % Continuation unitl node / boundary node is found
    maxiter=200;
    iii=1;
    if ~isempty(intersect(path,end_node))
        found=1;
    else
        found=0;
    end
    while found~=1 && iii<=maxiter
        for i=1:size(element,1)
            idxx(i)=length(intersect(connection_nodes,element{i,1}))==2;
        end
        idxx=find(idxx);
        idxx=setdiff(idxx,element_path(end));
        element_path=[element_path,idxx];
        
        if direction ==1
            temp=[element{idxx,1}, element{idxx,1}];
            temp_idx=find(path(end)==temp);
            temp_idx=temp_idx(2);
            path=[path, temp(temp_idx-1)];
            connection_nodes=[temp(temp_idx-1) temp(temp_idx-2)];
        elseif direction==2
            temp=[element{idxx,1}, element{idxx,1}];
            temp_idx=find(path(end)==temp);
            temp_idx=temp_idx(1);
            path=[path, temp(temp_idx+1)];
            connection_nodes=[temp(temp_idx+1) temp(temp_idx+2)];
        else
            
        end
        
        if ~isempty(intersect( path,end_node )) 
            found=1;
        elseif ~(length(unique( path))==length(path))
            found=1;
            element_path=0;
            path=0;
            warning('At least on continuation was unsuccesfull');
        end
        iii=iii+1;
    end
    if length(r)>1
        path_temp{k,1}=path;
        path=[start_node ];
        element_path_temp{k,1}=element_path;
    end
    
   
end
if length(r)>1
    path=path_temp;
    element_path=element_path_temp;
end 

end

