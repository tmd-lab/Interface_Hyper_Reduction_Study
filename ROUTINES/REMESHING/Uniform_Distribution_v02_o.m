function [NewElements_idx,sym_flag] = Uniform_Distribution_v02_o(Elements_idx,x,y,varargin)
% Uniform_Distribution_v02: Creation of Uniform reduced mesh
% INPUTS:
%   Elements_idx : cell containing elements of original mesh (list of node indexe)
%   x: Vector of x Coordinates of Nodes	
%   y: Vector of y Coordinates of Nodes	
%   varargin    : Number of Elements of reduced mesh in 
%                  Nx - No. of Elements in x direction
%                  Ny - No. of Elements in y direction
%                  Nh - No. of Elements in radial hole direction
% OUTPUTS:
%   NewElements_idx : cell array containing elements of reduced mesh (list
%   of node index of input node set)
%   sym_flag  : Flag = 1 if interface mesh is not symetric


%% Determine number of elments origianl mesh
[corner_idx,bolt_idx,~]=CreateBRBInterface(x,y);

[nl1, ~]=Node2NodePath('path',corner_idx(1),corner_idx(2),Elements_idx,x,y);
[nl2, ~]=Node2NodePath('path',corner_idx(2),corner_idx(3),Elements_idx,x,y);
[nl3, ~]=Node2NodePath('path',corner_idx(3),corner_idx(4),Elements_idx,x,y);
[nl4, ~]=Node2NodePath('path',corner_idx(4),corner_idx(1),Elements_idx,x,y);
if length(nl1)~=length(nl3) || length(nl2)~=length(nl4)
    warning('Something with your mesh is wrong')
end

if length(nl1)>length(nl2)
    Nx_original=length(nl1)-1;
    Ny_original=length(nl2)-1;
else
    Nx_original=length(nl2)-1;
    Ny_original=length(nl1)-1;
end

% Find all bifurcation points
NumberOfElemts=zeros(size(x));
for p=1:size(Elements_idx,1)
    NumberOfElemts(Elements_idx{p,1})=NumberOfElemts(Elements_idx{p,1})+1;
end
bifurcation=find(NumberOfElemts>4);

% boundary_nodes=unique([nl1,nl2,nl3,nl4,cell2mat(bolt_idx)',bifurcation']);
bifurcation_connection=cell(length(bifurcation),length(bifurcation));
bifurcation_connection_hole=cell(length(bifurcation),1);

intesect_bound=[];
for i=1:length(bifurcation)
    boundary_nodes=unique([nl1,nl2,nl3,nl4,cell2mat(bolt_idx)',...
        bifurcation(setdiff(1:length(bifurcation),i))']);
    [inter_temp, ~]=Node2NodePath('all',bifurcation(i),boundary_nodes,Elements_idx,x,y);
    for p=1:size(inter_temp,1)
        intesect_bound=[intesect_bound,inter_temp{p,1}(end)];
        idx_found=find(bifurcation==inter_temp{p,1}(end));
        if ~isempty(idx_found)
        bifurcation_connection{i,idx_found}=inter_temp{p,1};
        end
        idx_found=find(cell2mat(bolt_idx)'==inter_temp{p,1}(end));
        if ~isempty(idx_found)
        bifurcation_connection_hole{i,1}=inter_temp{p,1};
        end
    end
end
if isequal(size(bifurcation_connection_hole{1}),size(bifurcation_connection_hole))
error('Your Mesh around all bolt holes are different. To allow different mesh updata code');
% Possible impovement. Allow different Mesh for each hole
end
Nhole_original=size(bifurcation_connection_hole{1},2)-1;


[~,idx1,~]=intersect(nl1,intesect_bound);
[~,idx2,~]=intersect(nl2,intesect_bound);
[~,idx3,~]=intersect(nl3,intesect_bound);
[~,idx4,~]=intersect(nl4,intesect_bound);

if ~all(idx1==idx3(end:-1:1)) || ~all(idx2==idx4(end:-1:1))
    warning('The edges are not symmetric - We will see what happens ;)');
end

% Determine Number of Elements for each section
if idx1(1)>idx1(end)
    nelements1=abs(diff([Ny_original+1 idx1' 1]));
else
    nelements1=abs(diff([1 idx1' Ny_original+1]));
end

if idx2(1)>idx2(end)
    nelements2=abs(diff([Nx_original+1 idx2' 1]));
else
    nelements2=abs(diff([1 idx2' Nx_original+1]));
end
%% Plot Original Mesh
figure ('Name', 'Original Mesh')
hold on
%     pause(10)
for ii=1:size(Elements_idx,1)
    Poly_Elements=polyshape(x(Elements_idx{ii,:}),y(Elements_idx{ii,:}));
    %     pause(0.01)
    plot(Poly_Elements)
end
axis equal
%% User imput for New Mesh numer of Elements
if ~isempty(varargin)
    Nx=varargin{1,1};
    Ny=varargin{1,2};
    Nhole=varargin{1,3};
else
    prompt = {['Enter number of elements in x-direction [',num2str(length(nelements2)),',',num2str(Nx_original),']:' ],...
        ['Enter number of elements in y-direction [',num2str(length(nelements1)),',',num2str(Ny_original),']:' ] ...
        ['Enter number of elements around hole in radial direction  [',num2str(1),',',num2str(Nhole_original),']:' ] };
    title = 'Input number of Elements for Mesh';
    dims = [1 80];
    answer = inputdlg(prompt,title); %,dims,definput)
    if ~isempty(answer)
        Nx=round(str2num(answer{1,1}));
        Ny=round(str2num(answer{2,1}));
        Nhole=round(str2num(answer{3,1}));
            
    else
        disp('No input')
        return
    end
end
if ~(Nx<0 || Nx>Nx_original || Ny<0 || Ny>Ny_original || Nhole>Nhole_original || Nhole<0)
    disp(['The following numbers were selected: Nx=',num2str(Nx),' , Ny=',num2str(Ny),' , Nhole=',num2str(Nhole)]);
else
    disp('Invalid input')
    return
end
%% Devide Interface in selected patchs (consideration of fixed nodes)
% Split interface in x Direction
num_split_x=ones(size(nelements2));
if Nx>length(nelements2)
    split1=floor(Nx/length(nelements2));
    modo1=mod(Nx,length(nelements2)); 
    if min(nelements2)<=split1
        num_split_x=ones(size(nelements2))*min(nelements2);
        modo1=modo1+((split1-min(nelements2))*length(nelements2));
    else
        num_split_x=ones(size(nelements2))*split1;
    end
    m=1;
    while m<=modo1
        max_val=max(nelements2./num_split_x);
        idx_max=find(max_val==(nelements2./num_split_x));
        if length(idx_max)==1
            num_split_x(idx_max)=num_split_x(idx_max)+1;
            m=m+1;
        elseif length(idx_max)<=modo1-m
            num_split_x(idx_max(1))=num_split_x(idx_max(1))+1;
            m=m+1;
        else % length(idx_max)>modo1-m
            num_dist=modo1-m+1; % Number of remaining elements to distribute
            if mod(num_dist,2)==0
                idx_mid=floor(length(idx_max)/2);
                for p=1:(num_dist/2)
                    num_split_x(idx_max(idx_mid-(p-1)))=num_split_x(idx_max(idx_mid-(p-1)))+1;
                    if mod(length(idx_max),2)==0   % length of idx_max is even
                        num_split_x(idx_max(idx_mid+(p-1)+1))=num_split_x(idx_max(idx_mid+(p-1)+1))+1;
                    else
                        num_split_x(idx_max(idx_mid+(p-1)+2))=num_split_x(idx_max(idx_mid+(p-1)+2))+1;
                    end
                end
            else
                idx_mid=round(length(idx_max)/2);
                num_split_x(idx_max(idx_mid))=num_split_x(idx_max(idx_mid))+1;
                for p=1:floor(num_dist/2)
                    factor=floor((length(idx_max)/(num_dist-1)));
                    num_split_x(idx_max(idx_mid-factor*p))=num_split_x(idx_max(idx_mid-factor*p))+1;
                    num_split_x(idx_max(idx_mid+factor*p))=num_split_x(idx_max(idx_mid+factor*p))+1;
                end
            end
            m=modo1+1;
        end
    end
end

% Split interface in y Direction
num_split_y=ones(size(nelements1));
if Ny>length(nelements1)
    split1=floor(Ny/length(nelements1));
    modo1=mod(Ny,length(nelements1));
    if min(nelements1)<=split1
        num_split_y=ones(size(nelements1))*min(nelements1);
        modo1=modo1+((split1-min(nelements1))*length(nelements1));
    else
        num_split_y=ones(size(nelements1))*split1;
    end
    m=1;
    while m<=modo1
        max_val=max(nelements1./num_split_y);
        idx_max=find(max_val==(nelements1./num_split_y));
        if length(idx_max)==1
            num_split_y(idx_max)=num_split_y(idx_max)+1;
            m=m+1;
        elseif length(idx_max)<=modo1-m
            num_split_y(idx_max(1))=num_split_y(idx_max(1))+1;
            m=m+1;
        else % length(idx_max)>modo1-m
            num_dist=modo1-m+1; % Number of remaining elements to distribute
            if mod(num_dist,2)==0
                idx_mid=floor(length(idx_max)/2);
                for p=1:(num_dist/2)
                    num_split_y(idx_max(idx_mid-(p-1)))=num_split_y(idx_max(idx_mid-(p-1)))+1;
                    if mod(length(idx_max),2)==0   % length of idx_max is even
                        num_split_y(idx_max(idx_mid+(p-1)+1))=num_split_y(idx_max(idx_mid+(p-1)+1))+1;
                    else
                        num_split_y(idx_max(idx_mid+(p-1)+2))=num_split_y(idx_max(idx_mid+(p-1)+2))+1;
                    end
                end
            else
                idx_mid=round(length(idx_max)/2);
                num_split_y(idx_max(idx_mid))=num_split_y(idx_max(idx_mid))+1;
                for p=1:floor(num_dist/2)
                    factor=floor((length(idx_max)/(num_dist-1)));
                    num_split_y(idx_max(idx_mid-factor*p))=num_split_y(idx_max(idx_mid-factor*p))+1;
                    num_split_y(idx_max(idx_mid+factor*p))=num_split_y(idx_max(idx_mid+factor*p))+1;
                end
            end
            m=modo1+1;
        end
    end
end

% Check if symetric mesh is possible
if  all(nelements1(1:round(length(nelements1)/2)) == nelements1(end:-1:round(length(nelements1)/2)))
    if  ~all(num_split_y(1:round(length(num_split_y)/2)) == num_split_y(end:-1:round(length(num_split_y)/2)))
        warning('Input values in y direction do not allow symmetric mesh. Try different input');
    end
end
if  all(nelements2(1:round(length(nelements2)/2)) == nelements2(end:-1:round(length(nelements2)/2)))
    if  ~all(num_split_x(1:round(length(num_split_x)/2)) == num_split_x(end:-1:round(length(num_split_x)/2)))
        warning('Input values in x direction do not allow symmetric mesh. Try different input');
    end
end


for sss=1:2
    if sss==1
        num_split=num_split_y;
        nelements=nelements1;
        idx=[1 ,idx3',Ny_original+1];
    else
        num_split=num_split_x;
        nelements=nelements2;
        idx=[1 ,idx2',Nx_original+1];
    end
    
    
    split=false(idx(end),1);
    
    for k=1:length(num_split)
        num=floor(nelements(k)/num_split(k));
        modo=mod(nelements(k),num_split(k));
%         idxadd=zeros(num_split(k)+1,1); % Initializeidxadd
        
        if modo==0
            idxadd=idx(k):num:idx(k+1);
        elseif modo~=0 && mod(nelements(k),2)==0  % Element Section to split is even
            if mod(modo,2)==0  % Modolo is even
                temp=round(linspace(idx(k),idx(k)+(idx(k+1)-idx(k))/2,num_split(k)/2+1));
                idxadd=[temp(1:end-1), idx(k+1)-fliplr(temp)+idx(k)];
            else  % Modolo is odd
                temp=round(linspace(idx(k),idx(k)+(idx(k+1)-idx(k))/2-1,num_split(k)/2+1));
                idxadd=[temp, idx(k+1)-fliplr(temp)+idx(k)];
            end
        elseif modo~=0 && mod(nelements(k),2)==1  % Element Section to split is odd
            if ~mod(modo,2)==0  % Modolo is odd
                % Determine side of interface
                idx_mid=round((idx(k+1)-idx(k))/2)+1;
                temp=round(linspace(1,idx_mid,num_split(k)/2+1));
                ttemp=round(linspace(idx_mid,nelements(k)+1,num_split(k)/2+1));
                idxadd=[temp(1:end-1),ttemp];
                
                if k>length(num_split)/2
                    idxadd=idx(k+1)-idxadd+1;
                elseif k<length(num_split)/2
                    idxadd=idxadd+idx(k)-1;
                else
                    error('This is unexpected :D');
                end
            else  % Modolo is even
                idxadd=round(linspace(idx(k),idx(k+1),num_split(k)+1)); 
            end
        end
        split(idxadd)=true;
    end
    %          y_split(idx(k)+idxadd-1)=true;
    if sss==1
        y_split=split;
    else
        x_split=split;
    end
end

% Check if split procedure was symetric (if possible)
sym_flag=0;
if  all(num_split_x(1:round(length(num_split_x)/2)) == num_split_x(end:-1:round(length(num_split_x)/2))) ...
        && all(num_split_y(1:round(length(num_split_y)/2)) == num_split_y(end:-1:round(length(num_split_y)/2)))
    if ~all(x_split(1:round(length(x_split)/2)) == x_split(end:-1:round(length(x_split)/2))) ...
            || ~all(y_split(1:round(length(y_split)/2)) == y_split(end:-1:round(length(y_split)/2)))
        sym_flag=1;
    end
end


% Split interface in radial bolt
split=round(linspace(1,Nhole_original+1,Nhole+1));
hole_split=false(Nhole_original+1,1);
hole_split(split)=true;


%% Create coarse Mesh - Outer boundary fisrt 
Temp_Elements_idx=Elements_idx;
NewElements_idx=cell(1,2);

nl1=nl1(y_split);
nl2=nl2(x_split);
nl3=nl3(y_split(end:-1:1));
nl4=nl4(x_split(end:-1:1));
nhole_original=cell2mat(bifurcation_connection_hole);
nhole=nhole_original(:,hole_split);

% Mesh refinement for first boundary side
nl1_temp=[];
for j=1:num_split_x(1)
    vertices_node=nl4(end-j);
    for i=1:length(nl1)-1
        nl1_temp(i)=vertices_node;
        [~,tempElem1]=Node2NodePath('path',nl1(i),vertices_node,Temp_Elements_idx,x,y);
        [~,tempElem2]=Node2NodePath('path',nl1(i),nl1(i+1),Temp_Elements_idx,x,y);
        
        [old_elem_temp,all_nodes,coner_temp]= BRB_Fill(Temp_Elements_idx,[tempElem1,tempElem2]);
        NewElements_idx{size(NewElements_idx,1)+1,1}=coner_temp;
        NewElements_idx{size(NewElements_idx,1),2}=all_nodes;
        for p=1:length(old_elem_temp)
            Temp_Elements_idx{old_elem_temp(p),1}(:)=[];
        end
        Temp_Elements_idx=Temp_Elements_idx(~cellfun('isempty',Temp_Elements_idx));
        vertices_node=setdiff(coner_temp,[vertices_node,nl1(i),nl1(i+1)]);
    end
    nl1_temp(i+1)=vertices_node;
    nl1=nl1_temp;
    nl1_temp=[];
end

% Mesh refinement for second boundary side
nl2_temp=[];
for j=1:num_split_y(end)
    vertices_node=nl1(end-j);
    for i=num_split_x(1)+1:length(nl2)-1
        nl2_temp(i)=vertices_node;
        [~,tempElem1]=Node2NodePath('path',nl2(i),vertices_node,Temp_Elements_idx,x,y);
        [~,tempElem2]=Node2NodePath('path',nl2(i),nl2(i+1),Temp_Elements_idx,x,y);
        
        [old_elem_temp,all_nodes,coner_temp]= BRB_Fill(Temp_Elements_idx,[tempElem1,tempElem2]);
        NewElements_idx{size(NewElements_idx,1)+1,1}=coner_temp;
        NewElements_idx{size(NewElements_idx,1),2}=all_nodes;
        for p=1:length(old_elem_temp)
            Temp_Elements_idx{old_elem_temp(p),1}(:)=[];
        end
        Temp_Elements_idx=Temp_Elements_idx(~cellfun('isempty',Temp_Elements_idx));
        vertices_node=setdiff(coner_temp,[vertices_node,nl2(i),nl2(i+1)]);
    end
    nl2_temp(i+1)=vertices_node;
    nl2=nl2_temp;
    nl2_temp=[];
end


% Mesh refinement for third boundary side
nl3_temp=[];
for j=1:num_split_x(end)
    vertices_node=nl2(end-j);
    for i=num_split_y(1)+1:length(nl3)-1
        nl3_temp(i)=vertices_node;
        [~,tempElem1]=Node2NodePath('path',nl3(i),vertices_node,Temp_Elements_idx,x,y);
        [~,tempElem2]=Node2NodePath('path',nl3(i),nl3(i+1),Temp_Elements_idx,x,y);
        
        [old_elem_temp,all_nodes,coner_temp]= BRB_Fill(Temp_Elements_idx,[tempElem1,tempElem2]);
        NewElements_idx{size(NewElements_idx,1)+1,1}=coner_temp;
        NewElements_idx{size(NewElements_idx,1),2}=all_nodes;
        for p=1:length(old_elem_temp)
            Temp_Elements_idx{old_elem_temp(p),1}(:)=[];
        end
        Temp_Elements_idx=Temp_Elements_idx(~cellfun('isempty',Temp_Elements_idx));
        vertices_node=setdiff(coner_temp,[vertices_node,nl3(i),nl3(i+1)]);
    end
    nl3_temp(i+1)=vertices_node;
    nl3=nl3_temp;
    nl3_temp=[];
end

% Mesh refinement for forth boundary side
nl4_temp=zeros(size(nl4));
for j=1:num_split_y(1)
    vertices_node=nl3(end-j);
    for i=num_split_x(end)+1:length(nl4)-1-num_split_x(1)
        nl4_temp(i)=vertices_node;
        [~,tempElem1]=Node2NodePath('path',nl4(i),vertices_node,Temp_Elements_idx,x,y);
        [~,tempElem2]=Node2NodePath('path',nl4(i),nl4(i+1),Temp_Elements_idx,x,y);
        
        [old_elem_temp,all_nodes,coner_temp]= BRB_Fill(Temp_Elements_idx,[tempElem1,tempElem2]);
        NewElements_idx{size(NewElements_idx,1)+1,1}=coner_temp;
        NewElements_idx{size(NewElements_idx,1),2}=all_nodes;
        for p=1:length(old_elem_temp)
            Temp_Elements_idx{old_elem_temp(p),1}(:)=[];
        end
        Temp_Elements_idx=Temp_Elements_idx(~cellfun('isempty',Temp_Elements_idx));
        vertices_node=setdiff(coner_temp,[vertices_node,nl4(i),nl4(i+1)]);
    end
    nl4_temp(i+1)=vertices_node;
    nl4=nl4_temp;
    nl4_temp=[];
end



%% Create coarse Mesh - between and around holes

%  Mesh refinement around an inbetween the holes
conected=zeros(length(bifurcation),1);
connection_completed=false(length(bifurcation),1);
idx_sort_con=false(length(bifurcation));
for kk=1:length(bifurcation)
    for ll=1:length(bifurcation)
        idx_sort_con(kk,ll)=~isempty(bifurcation_connection{kk,ll});
    end
end
[~,idx_sort_con]=sort(sum(idx_sort_con,2));


for k=1:length(bifurcation)
    kk=idx_sort_con(k);
    
    while connection_completed(kk)==0
        
        idx_con=false(length(bifurcation),1);
        for ll=1:length(bifurcation)
            idx_con(ll)=~isempty(bifurcation_connection{kk,ll});
        end
        idx_con=find(idx_con(:));
        
        
        if conected(kk)<2
            [~,tempElemhole]=Node2NodePath('path',bifurcation(kk),bifurcation_connection_hole{kk,1}(end),Temp_Elements_idx,x,y);
            for m=1:length(idx_con)
                if length(bifurcation_connection{kk,idx_con(m)})>1  % Reduce Computional Effort and avoid problem around holes for curved connections
                    [~,temptempElemcon]=Node2NodePath('path',bifurcation(kk),bifurcation_connection{kk,idx_con(m)}(2),Temp_Elements_idx,x,y);
                else
                    [~,temptempElemcon]=Node2NodePath('path',bifurcation(kk),bifurcation_connection{kk,idx_con(m)}(end),Temp_Elements_idx,x,y);
                end
                tempElemcon(m)=length(intersect(temptempElemcon,tempElemhole));
            end
            idx_connect=find(tempElemcon);
            idx_connect=idx_con(idx_connect(1));
            
            direction1_temp= bifurcation_connection{kk,idx_connect};
            direction2=nhole(kk,:);
            
            if (conected(kk)+ length(idx_con))<=2
                bifurcation_connection{kk,idx_connect}=[];
                conected(kk)=conected(kk)+1;
                bifurcation_connection{idx_connect,kk}=[];
                conected(idx_connect)=conected(idx_connect)+1;
            else
                conected(kk)=conected(kk)+1;
                conected(idx_connect)=conected(idx_connect)+1;
            end
            
            bolteconnection=true;
        elseif conected(kk)==2
            if length(idx_con)==2
                direction1_temp=bifurcation_connection{kk,idx_con(1)};
                direction2_temp=bifurcation_connection{kk,idx_con(2)};
                
                bifurcation_connection{kk,idx_con(1)}=[];
                bifurcation_connection{kk,idx_con(2)}=[];
                bifurcation_connection{idx_con(1),kk}=[];
                bifurcation_connection{idx_con(2),kk}=[];
                
                bolteconnection=false;
            else
                erorr('This is unexpected');
            end
            
        end

        int1=length(intersect( direction1_temp,nl1));
        int2=length(intersect( direction1_temp,nl2));
        int3=length(intersect( direction1_temp,nl3));
        int4=length(intersect( direction1_temp,nl4));
        
        if int1>=2
            idx_start=find(nl1==direction1_temp(1));
            idx_end=find(nl1==direction1_temp(end));
            if idx_start>idx_end
                direction1=nl1(idx_start:-1:idx_end);
            else
                direction1=nl1(idx_start:idx_end);
            end
        elseif int2>=2
            idx_start=find(nl2==direction1_temp(1));
            idx_end=find(nl2==direction1_temp(end));
            if idx_start>idx_end
                direction1=nl2(idx_start:-1:idx_end);
            else
                direction1=nl2(idx_start:idx_end);
            end
        elseif int3>=2
            idx_start=find(nl3==direction1_temp(1));
            idx_end=find(nl3==direction1_temp(end));
            if idx_start>idx_end
                direction1=nl3(idx_start:-1:idx_end);
            else
                direction1=nl3(idx_start:idx_end);
            end
        elseif int4>=2
            idx_start=find(nl4==direction1_temp(1));
            idx_end=find(nl4==direction1_temp(end));
            if idx_start>idx_end
                direction1=nl4(idx_start:-1:idx_end);
            else
                direction1=nl4(idx_start:idx_end);
            end
        elseif int1==1 && int3==1
            error('This case should not occur for our problem. If it does implment somethin smart')
        elseif int2==1 && int4==1
            if ~isempty(intersect(direction1_temp(1),nl2))
                direction1=direction1_temp(y_split(idx1(1):-1:idx1(2)));
            else
                direction1=direction1_temp(y_split(idx1(2):idx1(1)));
            end
        else
            error('This case should not occur for our problem. If it does implment somethin smart')
        end
                
        if  bolteconnection==0
            
            int1=length(intersect( direction2_temp,nl1));
            int2=length(intersect( direction2_temp,nl2));
            int3=length(intersect( direction2_temp,nl3));
            int4=length(intersect( direction2_temp,nl4));
            
            if int1>=2
                idx_start=find(nl1==direction2_temp(1));
                idx_end=find(nl1==direction2_temp(end));
                if idx_start>idx_end
                    direction2=nl1(idx_start:-1:idx_end);
                else
                    direction2=nl1(idx_start:idx_end);
                end
            elseif int2>=2
                idx_start=find(nl2==direction2_temp(1));
                idx_end=find(nl2==direction2_temp(end));
                if idx_start>idx_end
                    direction2=nl2(idx_start:-1:idx_end);
                else
                    direction2=nl2(idx_start:idx_end);
                end
            elseif int3>=2
                idx_start=find(nl3==direction2_temp(1));
                idx_end=find(nl3==direction2_temp(end));
                if idx_start>idx_end
                    direction2=nl3(idx_start:-1:idx_end);
                else
                    direction2=nl3(idx_start:idx_end);
                end
            elseif int4>=2
                idx_start=find(nl4==direction2_temp(1));
                idx_end=find(nl4==direction2_temp(end));
                if idx_start>idx_end
                    direction2=nl4(idx_start:-1:idx_end);
                else
                    direction2=nl4(idx_start:idx_end);
                end
            elseif int1==1 && int3==1
                error('This case should not occur for our problem. If it does implment somethin smart')
            elseif int2==1 && int4==1
                if ~isempty(intersect(direction2_temp(1),nl2))
                    direction2=direction2_temp(y_split(idx1(1):-1:idx1(2)));
                else
                    direction2=direction2_temp(y_split(idx1(2):idx1(1)));
                end
            else
                error('This case should not occur for our problem. If it does implment somethin smart')
            end
        end
        
        % Combination
        nldirection2_temp=zeros(size(direction1));
        for j=2:size(direction2,2)
            vertices_node=direction2(j);
            for i=1:length(direction1)-1
                nldirection2_temp(i)=vertices_node;
                if size(direction2,2)>2 && length(direction1)>2  % Prevent Problems for low number of elements
                    [~,tempElem1]=Node2NodePath('path',direction1(i),vertices_node,Temp_Elements_idx,x,y);
                    [~,tempElem2]=Node2NodePath('path',direction1(i),direction1(i+1),Temp_Elements_idx,x,y);
                else
                    if size(direction2,2)>2 || length(intersect(direction2,bifurcation))==1
                        [~,tempElem1]=Node2NodePath('path',direction1(i),vertices_node,Temp_Elements_idx,x,y);
                    else
                        [~,ttempElem1]=Node2NodePath('path',direction1(i),direction2_temp(round(length(direction2_temp)/2)),Temp_Elements_idx,x,y);
                        [~,tttempElem1]=Node2NodePath('path',direction2_temp(round(length(direction2_temp)/2)),vertices_node,Temp_Elements_idx,x,y);
                        tempElem1=[ttempElem1,tttempElem1];
                    end
                    if length(direction1)>2
                        [~,tempElem2]=Node2NodePath('path',direction1(i),direction1(i+1),Temp_Elements_idx,x,y);
                    else
                        [~,ttempElem2]=Node2NodePath('path',direction1(i),direction1_temp(round(length(direction1_temp)/2)),Temp_Elements_idx,x,y);
                        [~,tttempElem2]=Node2NodePath('path',direction1_temp(round(length(direction1_temp)/2)),direction1(i+1),Temp_Elements_idx,x,y);
                        tempElem2=[ttempElem2,tttempElem2];
                    end
                    if length(direction1)==2
                        [~,~,coner_temp]= BRB_Fill(Temp_Elements_idx,[tempElem1,ttempElem2]);
                        new_direction1_temp=setdiff(coner_temp,[direction1_temp(round(length(direction1_temp)/2)),direction1(i),vertices_node]);
                        direction1_temp=zeros(3,1);
                        direction1_temp(2)=new_direction1_temp;
                    elseif size(direction2,2)==2 && length(intersect(direction2,bifurcation))~=1
                        [~,~,coner_temp]= BRB_Fill(Temp_Elements_idx,[ttempElem1,tempElem2]);
                        new_direction2_temp=setdiff(coner_temp,[direction2_temp(round(length(direction2_temp)/2)),direction1(i),direction1(i+1)]);
                        direction2_temp=zeros(3,1);
                        direction2_temp(2)=new_direction2_temp;
                    end
                end
                
                [old_elem_temp,all_nodes,coner_temp]= BRB_Fill(Temp_Elements_idx,[tempElem1,tempElem2]);
                NewElements_idx{size(NewElements_idx,1)+1,1}=coner_temp;
                NewElements_idx{size(NewElements_idx,1),2}=all_nodes;
                for p=1:length(old_elem_temp)
                    Temp_Elements_idx{old_elem_temp(p),1}(:)=[];
                end
                Temp_Elements_idx=Temp_Elements_idx(~cellfun('isempty',Temp_Elements_idx));
                vertices_node=setdiff(coner_temp,[vertices_node,direction1(i),direction1(i+1)]);
            end
            nldirection2_temp(i+1)=vertices_node;
            direction1=nldirection2_temp;
            nldirection2_temp=[];
        end
        
        % Finde last corner node for area in between the bolts
        if  bolteconnection==0
            mm=find(vertices_node==bifurcation);
            bifurcation_connection{mm,idx_con(1)}=[];
            bifurcation_connection{idx_con(1),mm}=[];
        end
        % Determine which bifurcation points are fully conected
        for mm=1:length(bifurcation)
            for pp=1:length(bifurcation)
                temp(mm,pp)=isempty(bifurcation_connection{mm,pp});
            end
        end
        if any(all(temp)==1)
            connection_completed(all(temp))=1;
        end
    end
end
%% Clean New Element Vector up and creat lost of convex hull
NewElements_idx(1,:)=[];
for ii=1:size(NewElements_idx,1)
    idx_convex=convhull(x(NewElements_idx{ii,1}),y(NewElements_idx{ii,1}));
idx_convex=idx_convex(1:length(NewElements_idx{ii,1}));
NewElements_idx{ii,1}=NewElements_idx{ii,1}(idx_convex);
end

%% Determine ID of interior nodes for each new element
toleranceCoincident=1e-9;
[NewElements_idx] = FindInteriorNodes(NewElements_idx,[x,y],[x,y],toleranceCoincident);

disp(['The new created mesh consits of ',num2str(size(NewElements_idx,1)),' Elements']);

end

