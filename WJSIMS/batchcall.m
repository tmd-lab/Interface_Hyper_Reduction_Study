clc
clear all

sel_method = 'U';
tot = 15;

% sel_method = 'P';
% tot = 7;

% sel_method = 'PD';
% tot = 7;

for i=2:tot
    batch(@transient_wjrom_call, 0, {sel_method, i});
end
