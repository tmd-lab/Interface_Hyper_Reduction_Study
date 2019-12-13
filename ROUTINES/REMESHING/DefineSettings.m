function [Setting] = DefineSettings
% Settings for Pressure Based Meshing algorithm are defined 
%% Limits for Loop optimization
Setting.scale_min=.01;           % Scale Factor for Density of nodes on contours lower and upper limit
Setting.scale_max=1;
Setting.maxiter=10;              % Maximal number of iterations to find correct number of Elements

%% Parameter for Node placement
% Number of Elements for Surface fit, which is used to create contours
Setting.NumberofElementsContour=400;
% Contours withs length below with a factor below a this value comparted to
% maximal length are removed
Setting.RemoveContourLengthFactor=0.01;

%% Parameter for Node Postprocessing 
% Tolerance for Move Node 2 boundary ( Factor of Node Densitiy)
Setting.MoveBoundaryToleranceOut=1/8;
Setting.MoveBoundaryToleranceHoles=1/8;
% Tolerance to remove coincident placed nodes ( Factor of Node Densitiy)
Setting.RemoveCoincident=1/20;

% Parameter Equal Density Function
Setting.ED.maxiter=25;                  % Maximal number of iterations (No ill itteration for bad inpur parameters)
Setting.ED.checkArea=true;              % Add nodes based on Area
Setting.ED.checkAspectRatio=false;      % Add nodes based on Aspect Ratio
Setting.ED.checkEdgeLengthRatio=true;   % Add nodes based on Edge Length Ratio
Setting.ED.checkEquiangularSkew=false;   % Add nodes based on Equiangular Skewnes of Triangulation
Setting.ED.maxAreaFactor=5;             % Factor which has to be exceeded for node placement (relative to mean)
Setting.ED.maxAspectRatio=3.5;          % Factor which has to be exceeded for node placement (relative to mean)
Setting.ED.edgeLengthRatio=2.5;         % Factor which has to be exceeded for node placement (relative to mean)
Setting.ED.EquiangularSkew=0.75;         % Value which has to be exceeded for node placement (absolute criteria)

%% Settings for Mesh
 Setting.MeshMethod=2;                 % 1: Aspect Ratio; 2: Skewness; 3: Product of AR & Skewness
 Setting.MeshMethodLimit=Setting.ED.EquiangularSkew;          % Limit Mesh method ( limit=0: triangulation) (Elements with quality below limit stay as created by triangulation)

end

