function [ TopInterface, BottomInterface , Element, Input ] = ImportNodes (sel_preload,sel_torque,PathName)
% ImportNodes: imports interface nodes and load values for both contact
% surfaces. If possible the elmenets are also imported
% INPUTS:
%   sel_preload : Pressure or Force load are imported ['Force','Pressure']
%   sel_torque	: Select torque level ['7NM' , '20NM']
%   PathName    : path were the import files are locaed
% OUTPUTS:
%   TopInterface    : Struct with ID,x-,y-,z-Coordinate,preload for top interface NnxNn Adjacency Matrix
%   BottomInterface : Struct with ID,x-,y-,z-Coordinate,preload for bottom interface NnxNn Adjacency Matrix
%   Element         : Struct containing Elements from top and Bottom
%   Input           : Raw Data table read input
if strcmp(sel_preload,'Force')
    FileName=[['TOPCOORDFORCE',sel_torque,'.dat'];['BOTCOORDFORCE',sel_torque,'.dat']];  %File for top first, bottom second
    disp('Nodal Force values are imported');
elseif strcmp(sel_preload,'Pressure')
    FileName=[['TOPCOORDCPRESS',sel_torque,'.dat'];['BOTCOORDCPRESS',sel_torque,'.dat']];  %File for top first, bottom second
    disp('Nodal Pressure values are imported');
else
    error('implement');
end


%  Make sure folder was selected
if PathName==0
    disp('No folder was selected. Script was stopped!')
    return
end
Input.Top=table2array(readtable([PathName,FileName(1,:)],'Delimiter',' '));
Input.Bottom=table2array(readtable([PathName,FileName(2,:)],'Delimiter',' '));

TopInterface.ID=Input.Top(:,1);
TopInterface.x=Input.Top(:,2);
TopInterface.y=Input.Top(:,3);
TopInterface.z=Input.Top(:,4);
TopInterface.preload=Input.Top(:,5);

if all(TopInterface.preload<=0)
    TopInterface.preload=-TopInterface.preload;
end

BottomInterface.ID=Input.Bottom(:,1);
BottomInterface.x=Input.Bottom(:,2);
BottomInterface.y=Input.Bottom(:,3);
BottomInterface.z=Input.Bottom(:,4);
BottomInterface.preload=Input.Bottom(:,5);
if all(BottomInterface.preload<=0)
    BottomInterface.preload=-BottomInterface.preload;
end

if exist ([PathName,'IntElements.dat'],'file')==2
    Input.Element=table2array(readtable([PathName,'IntElements.dat'],'Delimiter',' '));
    Element.ID=Input.Element(:,1);
    Element.TopID=Input.Element(:,2:5);
    Element.BottomID=Input.Element(:,6:9);
    Element.Row=Input.Element(:,10:13);
    
    % Check if imported data are consistent
    if ~all(all(Element.TopID==TopInterface.ID(Element.Row)))
        warning('Top interface IDs and Element rows are not consistent');
    elseif ~all(all(Element.BottomID==BottomInterface.ID(Element.Row)))
        warning('Bottom interface IDs and Element rows are not consistent');
    end
else
    Element.ID=0;
end

if ~isempty(intersect(TopInterface.ID,BottomInterface.ID))
    warning('At least one Node is part of both contact surfaces.');
elseif  length(unique(TopInterface.ID)) ~= length(TopInterface.ID)
    warning('TopInterface Node list is not unique');
elseif length(unique(BottomInterface.ID)) ~= length(BottomInterface.ID)
    warning('BottomInterface Node list is not unique');
end

end

