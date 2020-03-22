function [corners_idx,bolt_idx,pgon] = CreateBRBInterface(x,y)
% CreateBRBInterface: Determines corner holes, bolt holes
% INPUTS:
%   x : Vector containing all x - Coordinates
%   y : Vector containing all y - Coordinates
% OUTPUTS:
%   corners_idx     : Index of corners of BRB  
%   bolt_idx        : cell, Index of boundary nodes at bolt hole of BRB 
%   pgon            : Polygon object of the BRB
    Nbolts=3;
    % bolt_dist=30.1752*1e-3;  %[m]
    bolt_dist = range(x)/4;
    tolerance_hole_mesh=0.5*1e-3; %[m]

    rx=[min(x) max(x)];
    ry=[min(y) max(y)];

    bnd=[rx ry]; %data bounds
    corners=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]); %data boundary corners
    corners_idx=[];
    for i=1:4
        [~,idx]=min( sqrt((x-corners(i,1)).^2+(y-corners(i,2)).^2));
        corners_idx=[corners_idx;idx];
    end
    % Sort Boundary nodes 
    sort    =   convhull(x(corners_idx),y(corners_idx));
    sort    =   sort(1:end-1);
    corners_idx  =  corners_idx(sort,:);

    % Calculate center of Bolt holes and determine idex of nodes at hole
    dx=(rx(2)+rx(1))/2;
    dy=(ry(2)+ry(1))/2;
    for i=1:Nbolts
        bolt_center(i,:)=[corners(4,1)+dx/2+i*bolt_dist,dy];
        dist=sqrt((x-bolt_center(i,1)).^2+(y-bolt_center(i,2)).^2);
        [~,idx]=min(dist);
        bolt_idx{i,1}=find(dist<=dist(idx)+tolerance_hole_mesh);

        % Sort  nodes 
        sort    =   convhull(x(bolt_idx{i,1}(:)),y(bolt_idx{i,1}(:)));
        sort    =   sort(1:end-1);
        bolt_idx{i,1}  =  bolt_idx{i,1}(sort);
    end

    % Create Polygon with holes
    pgon=polyshape(x(corners_idx),y(corners_idx));
    for i=1:3
        pgon = addboundary(pgon,x(bolt_idx{i,1}(:)),y(bolt_idx{i,1}(:)));
    end

    % Check if all nodes are on polygon
    [in,on]=inpolygon(x,y,pgon.Vertices(:,1),pgon.Vertices(:,2)); % check if all points are in / on edge of polygon;
    if ~length(find((in + on)>0))==length(x)
        warning('Not all nodes are located on the interface');
    else
        %     disp('All nodes are located at the interface');
    end

    % figure
    % hold on
    % plot (x,y,'+');
    % plot(corners(:,1),corners(:,2),'r-');
    % for i=1:3
    % plot(x(bolt_idx{i,1}(:,1)),y(bolt_idx{i,1}(:,1)),'b-');
    % end
    % plot(pgon)
    % axis equal
    % axis tight
end

