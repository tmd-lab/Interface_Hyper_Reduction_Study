function [Elements_idx,qualityflag] = MeshGeometryBased(New_Nodes,methodIDX,selection_methodLimit,idx_holeBound)
% MeshGeometryBased generates mesh based on node input
% INPUTS:
%   New_Nodes: Array with Coordinates mesh	
%   methodIDX: meshing criteria 1:AspectRatio 2: Skew 3: AS & Skew 
%   selection_methodLimit: Limit meshin criteria - Elements with higher standard are generated (if possible)
%              Vector with preload (eg. Force / Pressure)
%   idx_holeBound : index of all nodes that are part of boundary along the hole          NumberOfElements - "Target" No. of Elements for reduced
% OUTPUTS:
%   NewElements_idx : cell array containing elements of reduced mesh (list
%   of node index of NewNodelList), and interior nodes of original mesh for
%   each element
%   qualityflag: If meshing criteria can not be fulfilled -->1

%% Create triangle and quadrilateral mesh
qualityflag=0;

% Create Mesh
selection_method=string([string('AspectRatio'),string('Skew'),string('AS_Skew')]);
selection_method=selection_method(methodIDX);
maxiter=25; % Just needed as backup - normally never needed

DT=delaunayTriangulation(New_Nodes);  % Delaunay Triangulation
Connectivity_Red=DT.ConnectivityList;
NTriangles = size(DT.ConnectivityList,1);

%  Remove Elements in Bolt holes
delete=false(NTriangles,1);
for i=1:length(idx_holeBound)
    delete=or(delete,DT.ConnectivityList==idx_holeBound(i));
end
Connectivity_Red(all(delete,2),:)=[];
NTri = size(Connectivity_Red,1);
% Creat triangulation with reduced connectivity - (Reason: DT is read only
% - connections in holes can not be removed
TR=triangulation(Connectivity_Red,New_Nodes); 
N = neighbors(TR);    % Find neighbor for each elemnt

% Establish strct which contains properties for possibly combined triangles
CEProperties=cell(NTri,3); 
EProperties=cell(NTri,1);
for k=1:NTri
     % Determine Geometric properties of Triangle
      EProperties{k}.Element.isBound=any(isnan(N(k,:)));
      EProperties{k}.Element.Nodes=Connectivity_Red(k,:);
      EProperties{k}.Element.Neighbors=N(k,:);
      % Edge Length
      coords=New_Nodes(Connectivity_Red(k,:),:);
      coords=diff([coords;coords(1,:)]);
      for l=1:length(coords)
          EProperties{k}.Element.EdgeLength(l)=norm(coords(l,:));
      end
      EProperties{k}.Element.AspectRatio=max(  EProperties{k}.Element.EdgeLength) /...
          min(  EProperties{k}.Element.EdgeLength);
      % Angle between sides -> Equiangular skew
      coords=New_Nodes(Connectivity_Red(k,:),:);
      coordAngle=diff([coords(end,:);coords;coords(1,:)]);
      for p=2:length(coordAngle)
          u= - coordAngle(p-1,:);
          v= coordAngle(p,:);
          EProperties{k}.Element.Theta(p-1)= acosd(dot(u,v)/(norm(u)*norm(v)));
      end
      ThetaMax=max(EProperties{k}.Element.Theta);
      ThetaMin=min(EProperties{k}.Element.Theta);
      ThetaE=60;   % for triangular Elements
      EProperties{k}.Element.EquiangularSkew= ...
          max([(ThetaMax-ThetaE)/(180-ThetaE),(ThetaE-ThetaMin)/(ThetaE)]);
      
      
      % Determine Geometric properties of Combination with each boundary
      for i=1:3
          CEProperties{k,i}=struct('Side',[]);
          if ~isempty(CEProperties{k,i}) % See if data are already avialable
              if isnan(N(k,i))
                  CEProperties{k,i}.Side.isBound=true;
              else
                  CEProperties{k,i}.Side.isBound=false;
                  CEProperties{k,i}.Side.ConElement=N(k,i);
                  % Nodes for combined Element
                  nodes=unique([Connectivity_Red(k,:), Connectivity_Red(N(k,i),:)]);
                  idx=convhull(New_Nodes(nodes,1),New_Nodes(nodes,2));
                  if length(idx)-1==length(nodes) % Test if new element is convex
                      CEProperties{k,i}.Side.isConvex=true;
                      CEProperties{k,i}.Side.Nodes= nodes(idx(1:end-1));
                      % Edge length for combined Element
                      coords=New_Nodes(nodes(idx),:);
                      coords=diff(coords);
                      for l=1:length(coords)
                          CEProperties{k,i}.Side.EdgeLength(l)=norm(coords(l,:));
                      end
                      CEProperties{k,i}.Side.AspectRatio=max( CEProperties{k,i}.Side.EdgeLength) /...
                          min( CEProperties{k,i}.Side.EdgeLength);
                      % Angle between sides -> Equiangular skew 
                      idxAngle=[idx(end-1);idx];
                      coordAngle=diff(New_Nodes(nodes(idxAngle),:));
                      for p=2:length(coordAngle)
                          u= - coordAngle(p-1,:);
                          v= coordAngle(p,:);
                          CEProperties{k,i}.Side.Theta(p-1)= acosd(dot(u,v)/(norm(u)*norm(v)));
                      end
                      ThetaMax=max(CEProperties{k,i}.Side.Theta);
                      ThetaMin=min(CEProperties{k,i}.Side.Theta);
                      ThetaE=90;   % for quadrilateral Elements
                      CEProperties{k,i}.Side.EquiangularSkew= ...
                          max([(ThetaMax-ThetaE)/(180-ThetaE),(ThetaE-ThetaMin)/(ThetaE)]);
                      else
                      CEProperties{k,i}.Side.isConvex=false;
                      CEProperties{k,i}.Side.Nodes=nodes;
                  end
                  
              end
              % Save Result for other Element too
              if isfield(CEProperties{k,i}.Side,'ConElement')
              [Idx]=find(N(CEProperties{k,i}.Side.ConElement,:)==k);
              CEProperties(CEProperties{k,i}.Side.ConElement,Idx)=CEProperties(k,i);
              CEProperties{CEProperties{k,i}.Side.ConElement,Idx}.Side.ConElement=k;

              end
          end
      end
end

%% Define Routine which selects order of "Best partner" - prevents connection

if strcmp(selection_method,'AspectRatio')
    Idx_sort=2;
elseif strcmp(selection_method,'Skew')
    Idx_sort=3;
elseif strcmp(selection_method,'AS_Skew')
    Idx_sort=4;
else
    error('Invalid Method selected');
end

% Determin ideal combination for each element 
for k=1:NTri
    Properties=zeros(4,3);
    del=[];
    for i=1:3
        if CEProperties{k,i}.Side.isBound==false && CEProperties{k,i}.Side.isConvex==true
            Properties(1,i)=round(CEProperties{k,i}.Side.ConElement);
            Properties(2,i)=abs(CEProperties{k,i}.Side.AspectRatio-1);
            Properties(3,i)=CEProperties{k,i}.Side.EquiangularSkew;
            Properties(4,i)=abs(CEProperties{k,i}.Side.AspectRatio-1)*CEProperties{k,i}.Side.EquiangularSkew;
        else
            del=[del;i];
        end
    end
%     Properties(1,4)=round(k);
%     Properties(2,4)=abs(EProperties{k}.Element.AspectRatio-1);
%     Properties(3,4)=EProperties{k}.Element.EquiangularSkew;
%     Properties(4,4)=abs(EProperties{k}.Element.AspectRatio-1)*EProperties{k}.Element.EquiangularSkew;
    
    Properties(:,del)=[];
    %Sort Matrix depending on selected parameter
    [~,idxSorted]=sort( Properties(Idx_sort,:),2);
    Properties=Properties(:,idxSorted);
    Properties=Properties(:,Properties(Idx_sort,:)<=selection_methodLimit);
    
    if isempty(Properties)
        
%     Add initial as Last possible Connection method, independently of geometry parameter
    Properties(1,1)=round(k);
    Properties(2,1)=abs(EProperties{k}.Element.AspectRatio-1);
    Properties(3,1)=EProperties{k}.Element.EquiangularSkew;
    Properties(4,1)=abs(EProperties{k}.Element.AspectRatio-1)*EProperties{k}.Element.EquiangularSkew;
    if EProperties{k}.Element.EquiangularSkew >selection_methodLimit
        qualityflag=1;
    end
    end
    EProperties{k}.Element.ConnectionPreference=Properties(1,:);
end

Elements_idx=cell(1,1);
isConnected=false(NTri,1);
Connection=zeros(NTri,2);
Elements_idx(1,:)=[];
kred=1:NTri;
p=1;
while all(isConnected)~= 1 && p<=maxiter
    % Connect Elements with same preference
    for k=1:length(kred)
        if EProperties{kred(k)}.Element.ConnectionPreference(1,1)~=kred(k) && isConnected(kred(k))==0
            idx_connection=EProperties{kred(k)}.Element.ConnectionPreference(1,1);
            if EProperties{idx_connection}.Element.ConnectionPreference(1,1)==kred(k)
                
                Connection(kred(k),1)=kred(k);
                Connection(kred(k),2)=idx_connection;
                Connection(idx_connection,1)=idx_connection;
                Connection(idx_connection,2)=kred(k);
                
                isConnected(kred(k))=true;
                isConnected(idx_connection)=true;
                
                % Add New Element to list of created Elements
                ConIdx=find(EProperties{kred(k)}.Element.Neighbors==idx_connection);
                NodeUnsorted=CEProperties{kred(k), ConIdx}.Side.Nodes;
                % Sort with respect to bigest x change
                [~,idx_Nsort]=max(diff([New_Nodes(NodeUnsorted,1);New_Nodes(NodeUnsorted(1),1)]));
                NodeUnsorted=[NodeUnsorted,NodeUnsorted];
                Elements_idx{end+1,1}=NodeUnsorted(idx_Nsort:idx_Nsort+3);

%                 Elements_idx{end+1,1}=CEProperties{kred(k), ConIdx}.Side.Nodes;
            end
            
        elseif EProperties{kred(k)}.Element.ConnectionPreference(1,1)==kred(k) && ...
                and( length(EProperties{kred(k)}.Element.ConnectionPreference(1,:))==1, isConnected(kred(k))==0)
            isConnected(kred(k))=true;
            Connection(kred(k),1)=kred(k);
            Connection(kred(k),2)=kred(k);
            
            % Sort Nodes
            NodeUnsorted=EProperties{kred(k)}.Element.Nodes;
            % Sort with respect to bigest x change
            [~,idx_Nsort]=max(diff([New_Nodes(NodeUnsorted,1);New_Nodes(NodeUnsorted(1),1)]));
            NodeUnsorted=[NodeUnsorted,NodeUnsorted];
            Elements_idx{end+1,1}=NodeUnsorted(idx_Nsort:idx_Nsort+2);
%             Elements_idx{end+1,1}=EProperties{kred(k)}.Element.Nodes;
            
        else
            %         error('This should not happen');
        end
    end
    
    % Reduce Connection Matrix of Elements
    kred=find(isConnected==0);
    connectedElements=unique(Connection);
    
    for k=1:length(kred)
        [allready_connectedElements,idx_remove]=intersect(EProperties{kred(k)}.Element.ConnectionPreference(1,:) ,connectedElements);
         EProperties{kred(k)}.Element.ConnectionPreference(:,idx_remove)=[];
        if isempty(EProperties{kred(k)}.Element.ConnectionPreference)
            
            EProperties{kred(k)}.Element.ConnectionPreference(1,1)=round(kred(k));
            if EProperties{kred(k)}.Element.EquiangularSkew >selection_methodLimit
                qualityflag=1;
            end
        end
    end
    p=p+1;
end
if p>maxiter
    error('Meshing wrong');
end
end