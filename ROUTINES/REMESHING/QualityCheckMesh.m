function [EquiangularSkew,AR,Area] = QualityCheckMesh(OutputNodes,OutputQuad,OutputTri,flag)
% QualityCheckMesh: Returns Quality paremeter for each element of mesh (Skewness, Aspect Ratio and Area)
% INPUTS:
%   OutputNodes,OutputQuad,OutputTri: Output of "CreateOutputFiles" Function 	
%   flag:  Flag if Values are display in command window
% OUTPUTS:
%   EquiangularSkew : Vector with Skewness for each element
%   AR              : Vector with Aspect Ratio for each element
%   Area            : Vector with Area for each element

%% Quad Elements
AreaQuad=zeros(size(OutputQuad,1),1);
EquiangularSkewQuad=AreaQuad;
ARQuad=AreaQuad;
EdgeLengthQuad=zeros(size(OutputQuad,1),4);
ThetaLengthQuad=zeros(size(OutputQuad,1),4);

for k=1:size(OutputQuad,1)
   coord=OutputNodes(OutputQuad(k,2:end),(2:3));
   Polyob=polyshape(coord(:,1),coord(:,2));
   % Calculate Area
   AreaQuad(k)=area(Polyob);
%    plot(Polyob)
   
   coords=[coord;coord(1,:)];
   coords=diff(coords);
   EdgeLengthQuad(k,:)=vecnorm(coords,2,2);
   % Calculate Aspect Ratio
   ARQuad(k)=max(EdgeLengthQuad(k,:))/min(EdgeLengthQuad(k,:));

   % Angle between sides -> Equiangular skew
   idxAngle=[4 1 2 3 4 1];
   coordAngle=diff(coord(idxAngle,:));
   for p=2:length(coordAngle)
       u= - coordAngle(p-1,:);
       v= coordAngle(p,:);
       ThetaLengthQuad(k,p-1)= acosd(dot(u,v)/(norm(u)*norm(v)));
   end
   ThetaMax=max(ThetaLengthQuad(k,:));
   ThetaMin=min(ThetaLengthQuad(k,:));
   ThetaE=90;   % for quadrilateral Elements
   EquiangularSkewQuad(k)= max([(ThetaMax-ThetaE)/(180-ThetaE),(ThetaE-ThetaMin)/(ThetaE)]);

end


%% Triangular Elements
AreaTri=zeros(size(OutputTri,1),1);
EquiangularSkewTri=AreaTri;
ARTri=AreaTri;
EdgeLengthTri=zeros(size(OutputTri,1),3);
ThetaLengthTri=zeros(size(OutputTri,1),3);

for k=1:size(OutputTri,1)
   coord=OutputNodes(OutputTri(k,2:end),(2:3));
   Polyob=polyshape(coord(:,1),coord(:,2));
   % Calculate Area
   AreaTri(k)=area(Polyob);
%    plot(Polyob)
   
   coords=[coord;coord(1,:)];
   coords=diff(coords);
   EdgeLengthTri(k,:)=vecnorm(coords,2,2);
   % Calculate Aspect Ratio
   ARTri(k)=max(EdgeLengthTri(k,:))/min(EdgeLengthTri(k,:));

   % Angle between sides -> Equiangular skew
   idxAngle=[3 1 2 3 1];
   coordAngle=diff(coord(idxAngle,:));
   for p=2:length(coordAngle)
       u= - coordAngle(p-1,:);
       v= coordAngle(p,:);
       ThetaLengthTri(k,p-1)= acosd(dot(u,v)/(norm(u)*norm(v)));
   end
   ThetaMax=max(ThetaLengthTri(k,:));
   ThetaMin=min(ThetaLengthTri(k,:));
   ThetaE=60;   % for triangular Elements
   EquiangularSkewTri(k)= max([(ThetaMax-ThetaE)/(180-ThetaE),(ThetaE-ThetaMin)/(ThetaE)]);

end
% 
EquiangularSkew=[EquiangularSkewQuad;EquiangularSkewTri];
AR=[ARQuad;ARTri];
Area=[AreaQuad;AreaTri];
if flag ==1
    disp(['Maximal Equiangular Skewness: ' num2str(max(EquiangularSkew))]);
    disp(['Minimal Equiangular Skewness: ' num2str(min(EquiangularSkew))]);
    disp(['Maximal Aspect Ratio: ' num2str(max(AR))]);
    disp(['Minimal Aspect Ratio: ' num2str(min(AR))]);
    disp(['Area Ratio (Amax/Amin): ' num2str((max(Area)/min(Area)))]);
end

end

