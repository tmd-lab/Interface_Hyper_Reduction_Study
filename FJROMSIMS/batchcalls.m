clc
clear all

% sel_method = 'U';
% tot = 7;

% sel_method = 'P';
% tot = 6;

sel_method = 'PD';
tot = 6;

for i=3:tot
    batch(@transient_fjrom_call, 0, {sel_method, i});
end


