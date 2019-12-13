function [ X1,X2 ] = Uniform_Distribution_v03(TopInterface, BottomInterface,Nx,Ny)
% This function creates a uniformly distributed mesh for the BRB interface
% for different Number in Elements in each axis.
% Input:
%    - TopInterface:     Struct with Nodal information of Top surfaces
%    (ID, Coordinates (x,y,z))
%    - BottomInterface:  Struct with Nodal information of Bottom surfaces
%    (ID, Coordinates (x,y,z))
%    - Nx: Number of Elements in beam axis
%    - Ny: Number of Elements orthogonal to beam axis

% Loop over both surfaces
for sel_interface=1:1
    warning('Add Top interface :D ');
    
    if  sel_interface==1        %     Bottom surfaces
        ID=BottomInterface.ID;
        x=BottomInterface.x;
        y=BottomInterface.y;
        z=BottomInterface.z;
    elseif sel_interface==2     %     Top surfaces
        ID=TopInterface.ID;
        x=TopInterface.x;
        y=TopInterface.y;
        z=TopInterface.z;
    else
        error('There only two surfaces. Check your selection');
    end
    if Nx<=Ny
        warning('The number of patches in x direction is higer or eaqual to the number in x direction.');
    end
    %% Create polygon for interface and bolt holes.
    
    [corners_idx,bolt_idx,pgon]=CreateBRBInterface(x,y);
    
    triangle_idx=delaunay(x,y);   % Delaunay triangulation
    
    figure
    triplot(triangle_idx,BottomInterface.x,BottomInterface.y)
    
    %% Celan up triangles in the bolt holes
    all_bolt_idx= cell2mat(bolt_idx);
    
    delete_idx=logical(false(size(triangle_idx)));   % Initialize Matrix
    for i=1:size(all_bolt_idx,1)
        delete_idx=delete_idx+(all_bolt_idx(i)==triangle_idx);
    end
    a=5;
    delete_triangle=all(delete_idx,2);
    triangle_idx(delete_triangle==1,:)=[];
    triplot(triangle_idx,BottomInterface.x,BottomInterface.y)
    
    %Create polyshape object for Triangles
    rr=1;
    for ii=1:size(triangle_idx,1)
        poly_triangular(ii,1)=polyshape(x(triangle_idx(ii,:)),y(triangle_idx(ii,:)));
        if ~issimplified( poly_triangular(ii,1))
            warning('polyshap object is not well formed');
        end
%         poly_triangular_area(ii,1)=area(poly_triangular(ii,1));
        [~,on]=inpolygon(x(triangle_idx(ii,:)),y(triangle_idx(ii,:)),...
            pgon.Vertices(:,1),pgon.Vertices(:,2));
        poly_triangular_boundary(ii,1)= (sum(on)>1 && sum(on)<=3);
        if (sum(on)>1 && sum(on)<=3)
        bound_node{rr,:}=triangle_idx(ii,on);
        rr=rr+1;
        end
    end
    
    %Find boundary Elemtents 
    bound=find(poly_triangular_boundary==1);
    % Find attached Elements
    del=[];
    elements_attached{size(bound,1),1}=[];
    for j=1:size(bound,1)
        [r1,~]=find(triangle_idx(bound(j),:)==triangle_idx(:,1));
        [r2,~]=find(triangle_idx(bound(j),:)==triangle_idx(:,2));
        [r3,~]=find(triangle_idx(bound(j),:)==triangle_idx(:,3));
        row=sort([r1;r2;r3]);
        idx_search=find(row==bound(j));
        row(idx_search)=[];
        elements_attached{j,1}=row(find(diff(row)<=1e-15));
        
        % Falsely found boundary elements
        if size(elements_attached{j,1},1)>2
            del=[del, j];
        end
    end
    
    % Remove incorrectly found boundary elements
    bound(del)=[];
    for p=1:size(del,2)
        elements_attached{del(p),1} =[];
        bound_node{del(p),:}=[];
    end
    elements_attached=elements_attached(~cellfun('isempty',elements_attached));
    bound_node=bound_node(~cellfun('isempty',bound_node));

        
    % Select right Element for mergeing.
    combine=zeros(size(bound));
    for j=1:size(bound,1)
        % Determine longest side
        if size(elements_attached{j,1},1)==1
            combine(j)= elements_attached{j,1};
        else
            node_idx=setdiff(triangle_idx(bound(j),:),bound_node{j,1});
            a1=[x(bound_node{j,1}(1))-x(bound_node{j,1}(2)) y(bound_node{j,1}(1))-y(bound_node{j,1}(2))];
            b1=[x(bound_node{j,1}(1))-x(node_idx) y(bound_node{j,1}(1))-y(node_idx)];
            alpha1 = acos( dot(a1,b1) / (norm(a1) * norm(b1)) );
            if alpha1> pi/2
                alpha1=alpha1-pi/2;
            end
            a2=[x(bound_node{j,1}(2))-x(bound_node{j,1}(1)) y(bound_node{j,1}(2))-y(bound_node{j,1}(1))];
            b2=[x(bound_node{j,1}(2))-x(node_idx) y(bound_node{j,1}(2))-y(node_idx)];
            alpha2 = acos( dot(a2,b2) / (norm(a2) * norm(b2)) );
            if alpha2> pi/2
                alpha2=alpha2-pi/2;
            end
            [~,idx]=min(abs([alpha1 alpha2]-pi/4));
            if idx==1 && all(size(idx)==[1 1])
                con_nodes=[bound_node{j,1}(1) node_idx];
            else
                con_nodes=[bound_node{j,1}(2) node_idx];
            end
            
            for i=1:size(elements_attached{j,1})
                temp=intersect(triangle_idx(elements_attached{j,1}(i),:),con_nodes);
                if size(temp,2)==2 
                combine(j)=elements_attached{j,1}(i);
                end
            end
            clear temp con_nodes idx node_idx
        end
    end
%% Combine plots
 pgon_reduced=pgon;
 if ~(length(combine)==length(unique(combine)))
     warning('Combine not unique')
 end
for ii=1:length(combine)
    
    pgon_union(ii)=union(poly_triangular(bound(ii)),poly_triangular(combine(ii)));
    pgon_reduced=subtract(pgon_reduced,poly_triangular(bound(ii)));
    pgon_reduced=subtract(pgon_reduced,poly_triangular(combine(ii)));
    
   hold on
   plot(poly_triangular(bound(ii)))
   pause(0.1)
   plot(poly_triangular(combine(ii)))
   pause(0.1)
   plot(pgon_union(ii))
   % Creat node list with pointe convex list ;) 
end
figure   
plot(pgon_reduced)



figure
hold on
% plot(TopInterface.x,TopInterface.y,'r*')
plot(x,y,'+')
text(x,y,string([num2str(x),num2str(y),num2str(ID)]));


    
    
    
    
end


figure
hold on
plot (x,y,'+');
plot(xTarget,yTarget,'*');
for j=1:length(patch_vertex_idx)
    plot(polyshape(x(patch_vertex_idx{j,1}(:)),y(patch_vertex_idx{j,1}(:))));
end
axis equal
axis tight

figure
plot(x,y,'+')
X2=0;
X1=0;
end

