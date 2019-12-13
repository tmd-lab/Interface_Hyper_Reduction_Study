function [pgon_result,idx_result,triangle] = RemoveTriangle(pgon,x,y,triangle_idx,varargin)
% RemoveTriangle: Combines triangels to on boundary to rectangular elements
% INPUTS:
%   pgon : Polygon object of the BRB
%   x    : Vector containing all x - Coordinates
%   y    : Vector containing all y - Coordinates
%   triangle_idx: index matrix from delauny triangulation
%   varargin: Results from previous iteration (element list)
% OUTPUTS:
%   pgon_results  : Polygon object of the reduced BRB (unmeshed BRB)
%   idx_result    : array containing all rectangular elements
%   triangle      : reduced triangulation array

%Create polyshape object for Triangles and determine boundary nodes
rr=1;
if isempty(varargin)
    for ii=1:size(triangle_idx,1)
        poly_triangular(ii,1)=polyshape(x(triangle_idx(ii,:)),y(triangle_idx(ii,:)));
        if ~issimplified( poly_triangular(ii,1))
            warning('polyshap object is not well formed');
        end
        [~,on]=inpolygon(x(triangle_idx(ii,:)),y(triangle_idx(ii,:)),...
            pgon.Vertices(:,1),pgon.Vertices(:,2));
        poly_triangular_boundary(ii,1)= (sum(on)>1 && sum(on)<=3);
        if (sum(on)>1 && sum(on)<=3)
            bound_node{rr,:}=triangle_idx(ii,on);
            rr=rr+1;
        end
    end
else
    bound_nodes_mesh=[];
    for jj=1:length(varargin{1,1})
        bound_nodes_mesh=[bound_nodes_mesh; varargin{1,1}{jj,:}(:) ];
    end
    bound_nodes_mesh=unique(bound_nodes_mesh);
    
    for ii=1:size(triangle_idx,1)
        poly_triangular(ii,1)=polyshape(x(triangle_idx(ii,:)),y(triangle_idx(ii,:)));
        if ~issimplified( poly_triangular(ii,1))
            warning('polyshap object is not well formed');
        end
        int=intersect(triangle_idx(ii,:),bound_nodes_mesh);
        poly_triangular_boundary(ii,1)= (length(int)>1 && length(int)<=3);
        if (length(int)>1 && length(int)<=3)
            bound_node{rr,:}=int;
            rr=rr+1;
        end
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
    hold on
    %     plot(poly_triangular(bound(j)))
    node_idx=setdiff(triangle_idx(bound(j),:),bound_node{j,1});
    % Determine longest side
    if size(elements_attached{j,1},1)==1
        combine(j)= elements_attached{j,1};
    elseif isempty(elements_attached{j,1})
        combine(j)=bound(j);
    elseif ~(length(node_idx)==1)
        num_att=zeros(size(elements_attached{j,1},1),1);
        for k=1:size(elements_attached{j,1},1)
            elem_att_idx=elements_attached{j,1}(k);
            m=1;
            stop=false;
            idx=setdiff(1:length(bound),j);
            while m<=size(bound,1)-1 && ~stop
                if ~isempty(find(elements_attached{idx(m),1}(:)==elem_att_idx))
                    num_att(k)=size(elements_attached{idx(m),1},1);
                    stop=true;
                end
                m=m+1;
            end
            
            if any(num_att==1)
                id=find(num_att==1);
                id=id(1);
                combine(j)=elements_attached{j,1}(id);
            elseif isempty(varargin)
                combine(j)=bound(j);
                warning('Check the find boundary nodes algorithm');
            else
                %%
                i=1;
                while i<=size(varargin{1,1},1) && ~(length(node_idx)==1)
                    if length(intersect(varargin{1,1}{i,:},bound_node{j,1}))==2
                        bound_node{j,1}= intersect(varargin{1,1}{i,:},bound_node{j,1});
                        node_idx=setdiff(triangle_idx(bound(j),:),bound_node{j,1});
                        %                    node_idx=setdiff(bound_node{j,1},intersect(varargin{1,1}{i,:},bound_node{j,1}));
                    end
                    i=i+1;
                end
            end
        end
    end
    if length(node_idx)==1
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
    end
    clear temp con_nodes idx node_idx a1 a2 b1 b2
end
%% Combine plots
if ~(length(bound)==length(unique(bound)))
    error('Something wrong with bound');
end

% Remove incongruent combinations
if ~(length([bound;combine])==length(unique([bound;combine])))
    
    [~,idx_int1,~]=intersect(bound,combine);
    [~,icombine,~]=unique(combine);
    idx_inc=setdiff(1:length(combine),icombine);
    idx_inc1=[];
    for k=1:length(idx_inc)
        temp1=find(combine==combine(idx_inc(k)));
        temp2=find(bound==combine(idx_inc(k)));
        idx_inc1=[idx_inc1;temp1 ; temp2];
    end
    index=[idx_int1;idx_inc1];
    
    %Check for inconsistency first
    
    del_idx=[];
    Elem_sort=sort([bound(index) combine(index)],2);
    [~,idx_sort]=unique(Elem_sort,'row');
    for k=1:length(idx_sort)
        elemtemp=Elem_sort(idx_sort(k),1);
        idx_temp=find(Elem_sort(:,1)==elemtemp);
        if ~(length(find(Elem_sort(:,2)==elemtemp))==0)
            if isempty(idx_temp)
                % Make sure single leftover elements are eliminated
                del=del;
            else
                del_idx=[del_idx ; unique([idx_temp;find(Elem_sort(:,2)==elemtemp)])];
            end
        elseif ~( length(unique(Elem_sort(idx_temp,2)))==1 )
            del_idx=[del_idx ; idx_temp];
        else
            for l=1:length(idx_temp)
                attach(l)=length(elements_attached{index(idx_temp(l)),1});
            end
            if min(attach)==1
                % Do not delete
            else
                [~,temp]=min(attach);
                del_idx=[del_idx ; idx_temp( setdiff(1:length(idx_temp),temp))];
            end
        end
    end
    bound(index(idx_sort(:)))=[];
    combine(index(idx_sort(:)))=[];
    for k=1:size(index(idx_sort(:)),1)
        bound_node{index(idx_sort(k)),1}=[];
        elements_attached{index(idx_sort(k)),1}=[];
    end
    bound_node=bound_node(~cellfun('isempty',bound_node));
    elements_attached=elements_attached(~cellfun('isempty',elements_attached));
end

if ~(length(combine)==length(unique(combine)))
    warning('Combine not unique')
elseif ~(length(union(combine,bound))==length(unique(union(combine,bound))))
    warning('Boundary cells and combine are not unique')
end
% figure(2)
% hold on

pgon_reduced=pgon;
for ii=1:length(combine)
    %      pgon_union(ii)=union(poly_triangular(bound(ii)),poly_triangular(combine(ii)));
    pgon_reduced=subtract(pgon_reduced,poly_triangular(bound(ii)));
    pgon_reduced=subtract(pgon_reduced,poly_triangular(combine(ii)));
    
    square_idx{ii,:}=unique([triangle_idx(bound(ii),:),triangle_idx(combine(ii),:)]);
    idx=convhull(x(square_idx{ii,1}),y(square_idx{ii,1}));
    idx=idx([1:end-1]);
    square_idx{ii,:}=square_idx{ii,:}(idx);
    
end
idx=union(bound,combine);
triangle_idx(idx,:)=[];

if isempty(varargin)
    Result=square_idx;
else
    Result=[varargin{1,1};square_idx];
end



pgon_result=pgon_reduced;
idx_result=Result;
triangle=triangle_idx;
end