function [New_Nodes] = EqualDensitiyNodes(Nodes,Setting)
% EqualDensitiyNodes: Adds additional nodes based on settings
% INPUTS:
%   Nodes:   List of nodes
%   Setting : Settings for addional node placement (e.g. criteria)
% OUTPUTS:
% New_Nodes: List of nodes with additionaly placed nodes

iteration=1;
endflag=0;
New_Nodes=Nodes;
AddNodes=[];

% Settings
maxiter=Setting.ED.maxiter; 
checkArea=Setting.ED.checkArea;
checkAspectRatio=Setting.ED.checkAspectRatio;
checkEdgeLengthRatio=Setting.ED.checkEdgeLengthRatio;
checkEquiangularSkew=Setting.ED.checkEquiangularSkew;
maxAreaFactor=Setting.ED.maxAreaFactor;
maxAspectRatio=Setting.ED.maxAspectRatio;
edgeLengthRatio=Setting.ED.edgeLengthRatio;
maxEquiangularSkew=Setting.ED.EquiangularSkew;

while endflag==0 && iteration<=maxiter
    DT=delaunayTriangulation(Nodes);
    NTriangles = size(DT.ConnectivityList,1);
% % %     Plot Triangulation
%     figure
%     triplot(DT)
%     axis equal
    
    % Determine Parameters for "Elements"
    Areas = zeros(NTriangles,1);
    EdgeLength=zeros(NTriangles,1);
    ElementSkewness=zeros(NTriangles,4);
    for i = 1:NTriangles
        PointIndexes = DT.ConnectivityList(i,:);
        % Area calculation
        Areas(i) = polyarea(DT.Points(PointIndexes,1),DT.Points(PointIndexes,2));
        polyin=polyshape(DT.Points(PointIndexes,1),DT.Points(PointIndexes,2));
        % Edge Length Caculation
        nodes_temp=[DT.Points(PointIndexes,:);DT.Points(PointIndexes(1),:)];
        nodes_temp=diff(nodes_temp);
        for k=1:3
            EdgeLength(i,k)=norm(nodes_temp(k,:));
        end
        %Skewness calculation
        coords=DT.Points(PointIndexes,:);
        coordAngle=diff([coords(end,:);coords;coords(1,:)]);
        for p=2:length(coordAngle)
            u= - coordAngle(p-1,:);
            v= coordAngle(p,:);
            ElementSkewness(i,p)= acosd(dot(u,v)/(norm(u)*norm(v)));
        end
        ThetaMax=max(ElementSkewness(i,2:4));
        ThetaMin=min(ElementSkewness(i,2:4));
        ThetaE=60;   % for triangular Elements
        ElementSkewness(i,1)= ...
            max([(ThetaMax-ThetaE)/(180-ThetaE),(ThetaE-ThetaMin)/(ThetaE)]);
        
      
    end
    MeanEdgeLength=mean(mean(EdgeLength));
    MeanArea=mean(Areas);
    
    % Add Nodes based on different criteria: Area, Aspect Ration, maximal
    % edge length
    AreaRatio=Areas/MeanArea;
    AspectRatio=max(EdgeLength,[],2)./min(EdgeLength,[],2);
    EdgeRatio=max(EdgeLength/MeanEdgeLength,[],2);
    if max(AreaRatio)>=maxAreaFactor  && checkArea==1 % Area based refinement
        [~,idx_maxArea]=max(Areas);
        ndirection= diff([DT.Points(DT.ConnectivityList(idx_maxArea,:),:);DT.Points(DT.ConnectivityList(idx_maxArea,1),:)]);
        AddNodes=DT.Points(DT.ConnectivityList(idx_maxArea,:),:)+ndirection/2;
        edgeRatio=EdgeLength(idx_maxArea,:)/MeanEdgeLength;
        AddNodes=AddNodes(edgeRatio==max(edgeRatio),:);
%         [~,idx_min]=min(edgeRatio);   % Alternative Definition
%         idx=setdiff(1:length(edgeRatio),idx_min);
%         AddNodes=mean(AddNodes(idx,:));
    elseif  max(EdgeRatio)>edgeLengthRatio && checkEdgeLengthRatio==1  % Certain edges are much longer than others
        [~,idx_maxER]=max(EdgeRatio);
        ndirection= diff([DT.Points(DT.ConnectivityList(idx_maxER,:),:);DT.Points(DT.ConnectivityList(idx_maxER,1),:)]);
        AddNodes=DT.Points(DT.ConnectivityList(idx_maxER,:),:)+ndirection/2;
        [~,idx_maxERi]=max(EdgeLength(idx_maxER,:),[],2);
        AddNodes=AddNodes(idx_maxERi,:);
    elseif max(ElementSkewness(:,1))>=maxEquiangularSkew  && checkEquiangularSkew==1 
        [~,idx_maxSkew]=max(ElementSkewness(:,1));
        ndirection= diff([DT.Points(DT.ConnectivityList(idx_maxSkew,:),:);DT.Points(DT.ConnectivityList(idx_maxSkew,1),:)]);
        AddNodes=DT.Points(DT.ConnectivityList(idx_maxSkew,:),:)+ndirection/2;
        edgeRatio=EdgeLength(idx_maxSkew,:)/MeanEdgeLength; % Add node in middel of longest side
        AddNodes=AddNodes(edgeRatio==max(edgeRatio),:);
    elseif max(AspectRatio) >maxAspectRatio && checkAspectRatio==1
        [~,idx_maxAR]=max(AspectRatio);
%         if AreaRatio(idx_maxAR)>=1 % Bad Aspect Ratio
            ndirection= diff([DT.Points(DT.ConnectivityList(idx_maxAR,:),:);DT.Points(DT.ConnectivityList(idx_maxAR,1),:)]);
            AddNodes=DT.Points(DT.ConnectivityList(idx_maxAR,:),:)+ndirection/2;
            [~,idx_maxARi]=max(EdgeLength(idx_maxAR,:),[],2);
            AddNodes=AddNodes(idx_maxARi,:);
%         else
%             AddNodes=[];
%         end
    else
        endflag=1;
    end
    iteration=iteration+1;
    if ~isempty(AddNodes)
        % Mirrow Nodes which are added
        
        
        New_Nodes=[New_Nodes; AddNodes];
    end
    Nodes=New_Nodes;
end
   
if iteration>maxiter+1
    warning('Node seed was ended by maximal number of iterations :(');
end

    
end

