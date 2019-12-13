function [NewElements_idx,NewNodelList] = PressureBasedDistribution(x,y,z,preload,varargin)
% PressureBasedDistribution: Creation of pressure /(preaload) based reduced mesh
% INPUTS:
%   x: Vector of x Coordinates of Nodes	
%   y: Vector of y Coordinates of Nodes	
%   z: Vector of z Coordinates of Nodes	
%   preload: Vector with preload (eg. Force / Pressure)
%   varargin    : Number of Elements of reduced mesh in 
%                  NumberOfElements - "Target" No. of Elements for reduced
%                  mesh
%                  ToleranceElements - Tolerance from "Target Value" 
% OUTPUTS:
%   NewElements_idx : cell array containing elements of reduced mesh (list
%   of node index of NewNodelList)
%   NewNodelList  : array with coordinates [x,y,z] of reduced nodes

  %% User imput for New Mesh numer of Elements and Tolerance
if ~isempty(varargin)
    NumberOfElements=varargin{1,1};
    ToleranceElements=varargin{1,2};
else
    prompt = {['Enter number of elements '],...
        ['Enter Tolerance for number of elements']};
        title = 'Input number of Elements for Mesh';
    dims = [1 80];
    answer = inputdlg(prompt,title); %,dims,definput)
    if ~isempty(answer)
        NumberOfElements=round(str2num(answer{1,1}));
        ToleranceElements=round(str2num(answer{2,1}));       
    else
        disp('No input')
        return
    end
end
 
if ToleranceElements<0
    error('Invalid input for Tolerance');
elseif  mod(NumberOfElements,4)>=ToleranceElements
%   ToleranceElements=4;
  warning('The choosen Tolerance is low - the algorithm may not converge');
end

% Input all Settings
[Setting] = DefineSettings;

scale.min=Setting.scale_min;
scale.max=Setting.scale_max;

%% Pressure based Mesh
[corners_idx,bolt_idx,~]=CreateBRBInterface(x,y);

% norm prelaod
preload=(preload-min(preload))/(max(preload)-min(preload));
% reduces model for to make sure result will be symmetric.
idx_reduced=intersect(find(x<=0),find(y<=0));
x_red=x(idx_reduced);
y_red=y(idx_reduced);
preload_red=preload(idx_reduced);

Nsegs = 5;  % 5 from Tobias
SegNumEle=zeros(1,Nsegs);
for N_Segments=1:Nsegs % Determine number of Segements which should be used
    
    scaleFactor=(scale.min +scale.max)/2;
    [New_Nodes,~,idx_holeBound,~] = PressureBasedNodePlacement_v02(x_red,y_red,x,y,preload_red,idx_reduced,bolt_idx,corners_idx,N_Segments,scaleFactor,Setting);
    [Elements_idx,qualityflag] = MeshGeometryBased(New_Nodes,Setting.MeshMethod, Setting.MeshMethodLimit,idx_holeBound);
    
    % % Number of Elements for based on different Number of segements for same scaling value
    SegNumEle(N_Segments)=length(Elements_idx)*4;
end
[~,N_Segments]=min(abs(SegNumEle-NumberOfElements));


% Loop to generate Mesh 

% Initial mesh generation with correct number of segemnts
scaleFactor=(scale.min +scale.max)/2;
[New_Nodes,~,idx_holeBound,dx] = PressureBasedNodePlacement_v02(x_red,y_red,x,y,preload_red,idx_reduced,bolt_idx,corners_idx,N_Segments,scaleFactor,Setting);
[Elements_idx,qualityflag] = MeshGeometryBased(New_Nodes,Setting.MeshMethod, Setting.MeshMethodLimit,idx_holeBound);
NEle =length(Elements_idx)*4;

% % plot
% figure;
% e3=find(cellfun(@(c) length(c), Elements_idx(:,1))==3);
% e4=find(cellfun(@(c) length(c), Elements_idx(:,1))==4);
% SHOW2DMESH(New_Nodes, [e3 cell2mat(Elements_idx(e3))], [e4 cell2mat(Elements_idx(e4))], 1, -1, -100); axis equal; axis off

iter=1;
maxiter_loop=Setting.maxiter;
while ~and(NEle<=NumberOfElements+ToleranceElements ,  NEle>=NumberOfElements-ToleranceElements) && iter<=  maxiter_loop
    if NEle-NumberOfElements>0
        scale.min=scaleFactor;
    else
        scale.max=scaleFactor;
    end
    scaleFactor=(scale.min +scale.max)/2;
    
    [New_Nodes,~,idx_holeBound,dx] = PressureBasedNodePlacement_v02(x_red,y_red,x,y,preload_red,idx_reduced,bolt_idx,corners_idx,N_Segments,scaleFactor,Setting);
    [Elements_idx,qualityflag] = MeshGeometryBased(New_Nodes,Setting.MeshMethod, Setting.MeshMethodLimit,idx_holeBound);
    
    % Multiplication, as only quarter Segement is analysed so far
    NEle=length(Elements_idx)*4;
    iter=iter+1;
end
% close all
if iter>maxiter_loop
warning('With the choosen settings no mesh could be generated. Check your input');
end
if qualityflag==1
     warning('At least one triagulation Element & all possible Connections do not fulfill your Mesh quality standard')
end

% Mirror in x Direction
BoundNodesidx=find(abs(New_Nodes(:,1)-dx)<=1e-8);
[New_Nodes,Elements_idx] = MirrorElements(New_Nodes,Elements_idx,BoundNodesidx,'X',dx);

% % plot
% figure;
% e3=find(cellfun(@(c) length(c), Elements_idx(:,1))==3);
% e4=find(cellfun(@(c) length(c), Elements_idx(:,1))==4);
% SHOW2DMESH(New_Nodes, [e3 cell2mat(Elements_idx(e3))], [e4 cell2mat(Elements_idx(e4))], 1, -1, -100); axis equal; axis off

% Mirror in y Direction
dy=0;
BoundNodesidx=find(abs(New_Nodes(:,2)-dy)<=1e-8);
[New_Nodes,Elements_idx] = MirrorElements(New_Nodes,Elements_idx,BoundNodesidx,'Y',dy);

% % plot
% figure;
% e3=find(cellfun(@(c) length(c), Elements_idx(:,1))==3);
% e4=find(cellfun(@(c) length(c), Elements_idx(:,1))==4);
% SHOW2DMESH(New_Nodes, [e3 cell2mat(Elements_idx(e3))], [e4 cell2mat(Elements_idx(e4))], 1, -1, -100); axis equal; axis off

%Determine minimal edge side
edgelengthmin=inf;
for k=1:length(Elements_idx)
    nodes=diff([New_Nodes(Elements_idx{k,1},1),New_Nodes(Elements_idx{k,1},2);New_Nodes(Elements_idx{k,1}(1),1),New_Nodes(Elements_idx{k,1}(1),2)]);
    edgelength=sqrt(nodes(:,1).^2+nodes(:,2).^2);
    edgelengthmin=min(edgelengthmin,min(edgelength));
end


%% Determine ID of interior nodes for each new element
toleranceCoincident=edgelengthmin/20;
[Elements_idx] = FindInteriorNodes(Elements_idx,New_Nodes,[x,y],toleranceCoincident);
NewElements_idx=Elements_idx;
% Add Coordinate
z_add=ones(length(New_Nodes),1)*z(1);
NewNodelList=[New_Nodes,z_add];
%%
disp(['The new created mesh consits of ',num2str(size(NewElements_idx,1)),' Elements']);
end




