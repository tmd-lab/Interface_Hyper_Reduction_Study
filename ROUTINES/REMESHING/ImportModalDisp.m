function [ ModXDispNorm ] = ImportModalDisp (ModNr,sel_torque,PathName)
% ImportModalDisp: imports displacements for spesific modeshape 
% surfaces. If possible the elmenets are also imported
% INPUTS:
%   ModNr       : Mode number
%   sel_torque	: Select torque level ['7NM' , '20NM']
%   PathName    : path were the import files are locaed
% OUTPUTS:
%   ModXDispNorm  : normalized modeshape 

% PathName = uigetdir('*.dat','Select folder with results of NL static analysis');
if ModNr==1
    FileName=['M1disps',sel_torque,'.dat'];  %File for first mode
elseif ModNr==2
    FileName=['M2disps',sel_torque,'.dat'];  %File for secon mode
else
    error('Invalid Input');
end

%  Make sure folder was selected
if PathName==0
    disp('No folder was selected. Script was stopped!')
    return
end
%% Import modal displacements


ModXDisp=table2array(readtable([PathName,FileName],'Delimiter',' '));
% Select DOF with Maximal displacement
MeanModXDisp=mean(ModXDisp);
ModXDisp=ModXDisp(:,abs(MeanModXDisp)==max(abs(MeanModXDisp)));
if all(ModXDisp<0)
ModXDisp=ModXDisp+abs(min(ModXDisp));
end

% Normalize preload and modal disp
ModXDispNorm=(ModXDisp-min(ModXDisp))/(max(ModXDisp)-min(ModXDisp));

end

