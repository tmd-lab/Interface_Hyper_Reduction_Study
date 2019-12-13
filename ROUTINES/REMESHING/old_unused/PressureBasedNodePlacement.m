function [New_Nodes,idx_OuterBound,idx_holeBound,dx] = PressureBasedNodePlacement(x_red,y_red,x,y,preload_red,idx_reduced,bolt_idx,corners_idx,N_Segments,scaleFactor,Setting)
% Determine new Boundary nodes for reduced interface
New_Nodes=[];
blc=[x(intersect(corners_idx,idx_reduced)),y(intersect(corners_idx,idx_reduced))]; % Bottom left corner
idx_temp=find(x_red==max(x_red));
[~,idx_min]=min(y_red(idx_temp));
idx_blc_red=idx_temp(idx_min);
brc=[x_red(idx_blc_red),y_red(idx_blc_red)]; % bottom right corner
idx_mid_hole=intersect(bolt_idx{2,1},idx_reduced);
[~,idx_min]=min(y(idx_mid_hole));
idx_mbhb=idx_mid_hole(idx_min);
mbhb = [x(idx_mbhb) ,y(idx_mbhb)];  % middle bolt hole bottom
[~,idx_max]=max(y(idx_mid_hole));
idx_mbht=idx_mid_hole(idx_max);
mbht=[x(idx_mbht) ,y(idx_mbht)];   % middle bolt hole top
idx_left_hole=intersect(bolt_idx{1,1},idx_reduced);
[~,idx_max]=max(x(idx_left_hole));
idx_lbhr=idx_left_hole(idx_max);
lbhr =[x(idx_lbhr) ,y(idx_lbhr)]; % left bolt hole right
[~,idx_min]=min(x(idx_left_hole));
idx_lbhl=idx_left_hole(idx_min);
lbhl =[x(idx_lbhl) ,y(idx_lbhl)]; % left bolt hole left
idx_temp=find(x_red==min(x_red));
[~,idx_max]=max(y_red(idx_temp));
idx_tlc_red=idx_temp(idx_max);
tlc=[x_red(idx_tlc_red),y_red(idx_tlc_red)]; % top left corner

clear idx_temp idx_min idx_max  idx_blc_red
clear idx_lbhr idx_lbhl idx_tlc_red idx_mbhb idx_mbht

[fit, ~] = CreateFit(x_red, y_red, preload_red);
split=linspace(min(preload_red),1,N_Segments+2);
split=split(2:end);
N_steps=Setting.NumberofElementsContour;
[x_mesh,y_mesh] = meshgrid(linspace(min(x_red),max(x_red),N_steps),linspace(min(y_red),max(y_red),N_steps));
preload_mesh=fit(x_mesh,y_mesh);
figure
[Contour_Lines,~] = contour(x_mesh,y_mesh,preload_mesh,split,'ShowText','on');
axis equal
% close

isoPreload=cell(1,3);
idx_start=1;
i=1;
while idx_start<= length(Contour_Lines)
    isoPreload{i,1}=Contour_Lines(1,idx_start);
    length_iso=Contour_Lines(2,idx_start);
    isoPreload{i,2}=Contour_Lines(2,idx_start);
    isoPreload{i,3}=Contour_Lines(:,idx_start+1:idx_start+length_iso);
    cont_length=0;
    for k=1:length_iso-1
        cont_length=cont_length+norm([(isoPreload{i,3}(1,k)-isoPreload{i,3}(1,k+1)) ...
            (isoPreload{i,3}(2,k)-isoPreload{i,3}(2,k+1))]);
    end
    isoPreload{i,4}=cont_length;
    idx_start=idx_start+length_iso+1;
    i=i+1;
end
% Ad holes to Contour

% Fit holes
% Middel hole
x_temp=x(idx_mid_hole(end))-cos(linspace(0,pi()/2,N_steps/2))*(x(idx_mid_hole(end))-x(idx_mid_hole(1))); % Equidistant placement of nodes
isoPreload{end+1,3}=[x_temp ;interp1(x(idx_mid_hole),y(idx_mid_hole),x_temp)];
isoPreload{end,3}(:,1)=mbht';
isoPreload{end,3}(:,end)=mbhb';
cont_length=0;
for k=1:length(isoPreload{end,3})-1
    cont_length=cont_length+norm([(isoPreload{end,3}(1,k)-isoPreload{end,3}(1,k+1)), ...
        (isoPreload{end,3}(2,k)-isoPreload{end,3}(2,k+1))]);
end
isoPreload{end,4}=cont_length;
isoPreload{end,2}=size(isoPreload{end,3},2);
%Left hole
x_temp=x(idx_left_hole(1))-(cos(linspace(0,pi(),N_steps/2))-1)*(x(idx_left_hole(end))-x(idx_left_hole(1)))/2; % Equidistant placement of nodes
isoPreload{end+1,3}=[x_temp ;interp1(x(idx_left_hole),y(idx_left_hole),x_temp)];
isoPreload{end,3}(:,1)=lbhl';
isoPreload{end,3}(:,end)=lbhr';
cont_length=0;
for k=1:length(isoPreload{end,3})-1
    cont_length=cont_length+norm([(isoPreload{end,3}(1,k)-isoPreload{end,3}(1,k+1)), ...
        (isoPreload{end,3}(2,k)-isoPreload{end,3}(2,k+1))]);
end
isoPreload{end,4}=cont_length;
isoPreload{end,2}=size(isoPreload{end,3},2);

% Points per length - will be iteration parameter (N_Segement, too)
ppl=0.0254/2*scaleFactor;

% Remove Contours which are too small compared to others
length_temp=cell2mat(isoPreload(:,4));
length_temp_norm=length_temp/(max(length_temp));
idx_del1=find(length_temp_norm(1:end-2)<=Setting.RemoveContourLengthFactor);
idx_del2=find(length_temp(1:end-2)/ppl<0.4);
idx_del=unique([idx_del1;idx_del2]);
if ~isempty(idx_del)
    isoPreload(idx_del,:)=[];
    length_temp(idx_del)=[];
    length_temp_norm(idx_del)=[];
end
    
    % Initialize boundarys to new node list
% New_Nodes=[blc;brc;mbhb;mbht;lbhr;lbhl;tlc];
New_Nodes=[blc;brc;tlc];

nNodes=length_temp/ppl+1;

if round(nNodes(end-1))<1 % Ensure Middle hole exists
    nNodes(end-1)=1;
end
if round(nNodes(end))<2 % Ensure left hole exists
    nNodes(end)=2;
end
% nNodes(end-1)=nNodes(end-1)+1;% add one node, as node will be dublicate for hole
nNodes(end)=nNodes(end)+.5;
%% Quarter Modell

for k=1:length(nNodes)
    idx=round(linspace(1,isoPreload{k,2},max(round(nNodes(k)+1),2)));
    New_Nodes=[New_Nodes;isoPreload{k,3}(:,idx)'];
end

%% Full model caclucaltion
% % % % for k=1:length(nNodes)
% % % %     % Determin on which side of contour vector a boundary is
% % % %     if isoPreload{k,3}(2,1)==blc(2) || isoPreload{k,3}(1,1)==blc(1)  % First node of vector is on lower boundary
% % % %         if isoPreload{k,3}(2,end) ==blc(2) || isoPreload{k,3}(1,end)==blc(1) % Both nodes are on outer boundary
% % % %             idx=round(linspace(1,isoPreload{k,2},max(round(nNodes(k)),2)));
% % % %             New_Nodes=[New_Nodes;isoPreload{k,3}(:,idx)'];
% % % %         else
% % % %             
% % % %             if mod(round(nNodes(k)*2),2)==0
% % % %                 idx=round(linspace(1,isoPreload{k,2},round(nNodes(k))+1));
% % % %             else
% % % %                 idx=round(linspace(1,2*isoPreload{k,2},round(nNodes(k)*2)+1));
% % % %                 iidx=find(idx<=isoPreload{k,2});
% % % %                 idx=idx(iidx);
% % % %             end
% % % %             New_Nodes=[New_Nodes;isoPreload{k,3}(:,idx)'];
% % % %         end
% % % %     elseif isoPreload{k,3}(2,end) ==blc(2) || isoPreload{k,3}(1,end)==blc(1) % End node is on outer boundary
% % % %         if mod(round(nNodes(k)*2),2)==0
% % % %             idx=round(linspace(isoPreload{k,2},1,round(nNodes(k))+1));
% % % %         else
% % % %             idx=round(linspace(1,2*isoPreload{k,2},round(nNodes(k)*2)+1));
% % % %             idx=isoPreload{k,2}-idx+1;
% % % %             iidx=find(idx>=1);
% % % %             idx=idx(iidx);
% % % %         end
% % % %         New_Nodes=[New_Nodes;isoPreload{k,3}(:,idx)'];
% % % %         % both nodes are on inner boundary
% % % %     elseif isoPreload{k,3}(2,1)==isoPreload{k,3}(2,end) || isoPreload{k,3}(1,1)==isoPreload{k,3}(1,end) %Both nodes are on same inner noundary - distribute
% % % % %         if mod(round(nNodes(k)*2),2)==1
% % % %             idx=round(linspace(1,isoPreload{k,2},round(nNodes(k))+1));
% % % % %         else
% % % % %             idx=round(linspace(1,2*isoPreload{k,2},round(nNodes(k)*2)+1));
% % % % %             iidx=find(idx<=isoPreload{k,2});
% % % % %             idx=idx(iidx);
% % % % %         end
% % % %         New_Nodes=[New_Nodes;isoPreload{k,3}(:,idx)'];
% % % %         
% % % %     elseif and( isoPreload{k,3}(1,1)==brc(1) , isoPreload{k,3}(2,end)==tlc(2) ) || ...
% % % %             and( isoPreload{k,3}(1,end)==brc(1) , isoPreload{k,3}(2,1)==tlc(2) ) % nodes are on different inner boundarys
% % % %         
% % % %         if mod(round(nNodes(k)*4),2)==0
% % % %             idx=round(linspace(1,2*isoPreload{k,2},round(nNodes(k)*2)+1));
% % % %             iidx=find(idx<=isoPreload{k,2});
% % % %             idx=idx(iidx);
% % % %         else %
% % % %             idx=round(linspace(1,2*isoPreload{k,2},round(nNodes(k)*2)+1));
% % % %             iidx=find(idx<=isoPreload{k,2});
% % % %             idx=idx(iidx);
% % % %         end
% % % % 
% % % % 
% % % %         New_Nodes=[New_Nodes;isoPreload{k,3}(:,idx)'];
% % % %     end
% % % % end

New_Nodes=unique(New_Nodes,'rows'); % Remove dublicates

New_Nodes = EqualDensitiyNodes(New_Nodes,Setting); % Add nodes for homogen mesh

% Create Polyshap object for hole areas - Remove nodes in hole  
polyh1=polyshape(x([idx_left_hole(1);idx_left_hole;idx_left_hole(end)]),[0.1;y(idx_left_hole);0.1]);
polyh2=polyshape([0.1;x(idx_mid_hole)],[0.1;y(idx_mid_hole)]);
[TFin1,TFon1] = isinterior(polyh1,New_Nodes);
[TFin2,TFon2] = isinterior(polyh2,New_Nodes);
TFin=or( xor(TFin1,TFon1) , xor(TFin2,TFon2));
% Remove nodes in holes
New_Nodes(TFin,:)=[];

% % Move nodes close to hole to boundary
tolerance=ppl*Setting.MoveBoundaryToleranceHoles;
idx_holeBound=[];
for k=1:size(bolt_idx,1)
    [New_Nodes,idx_temp] = Move2Boundary(New_Nodes,[x(bolt_idx{k,1}) y(bolt_idx{k,1})],tolerance);
    idx_holeBound=[idx_holeBound;idx_temp];
%     pgon = addboundary(pgon,New_Nodes(idx_temp,1),New_Nodes(idx_temp,2));
end



% % % % % Mirrow Nodes
dx=(x(corners_idx(1))+x(corners_idx(4)))/2; % due to offset in x-direction of (7.06e-4)
% % % % New_Nodes=[New_Nodes;(-(New_Nodes(:,1)-dx))+dx,New_Nodes(:,2)]; % x-Direction
% % % % New_Nodes=[New_Nodes;New_Nodes(:,1),-New_Nodes(:,2)]; % y-Direction
% % % % 
% % % % 
% % % % New_Nodes=unique(New_Nodes,'rows'); % Remove dublicates

%If a node is close to manully created boundary node remove this node
 [idx_remove_candidates, Distance]=rangesearch(New_Nodes,brc,ppl*Setting.RemoveCoincident*2);
 if length(idx_remove_candidates)==2
     New_Nodes(idx_remove_candidates(max(Distance)==Distance))=[];
 end
 
% Remove nodes, which almost coincident
tolerance=ppl*Setting.RemoveCoincident;
[ReducedNodes] = RemoveCloseNodes(New_Nodes,tolerance);
New_Nodes=ReducedNodes;



% Move nodes which are located almost at boundary
tolerance=ppl*Setting.MoveBoundaryToleranceOut;
% % For quarter Segment
[New_Nodes,idx_OuterBound] = Move2Boundary(New_Nodes,[tlc;blc;brc;brc(1),tlc(2)],tolerance);
% % For whole interface
% [New_Nodes,idx_OuterBound] = Move2Boundary(New_Nodes,[x(corners_idx),y(corners_idx)],tolerance);

% New_Nodes = EqualDensitiyNodes(New_Nodes); % Add nodes for homogen mesh

% % % % Find boundary nodes of holes
% % % [~,idx_holeBound,~] = CreateBRBInterface(New_Nodes(:,1),New_Nodes(:,2));
% % % idx_holeBound=cell2mat(idx_holeBound);

% Determine correct Idx for Bolt hole nodes
idx_holeBound=[];
for k=1:size(bolt_idx,1)
    [~,idx_temp] = Move2Boundary(New_Nodes,[x(bolt_idx{k,1}),y(bolt_idx{k,1})],1e-5);
    idx_holeBound=[idx_holeBound;idx_temp];
end

% % % New Nodes with ID
% figure
% hold on
% plot(New_Nodes(:,1),New_Nodes(:,2),'r*')
% Row=[1:length(New_Nodes(:,1))]';
% text(New_Nodes(:,1),New_Nodes(:,2),string([num2str(Row)]));
% axis equal

end

