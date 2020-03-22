clc
clear all

sel_method = 'U';
tot = 7;

for i=1:tot
    pp(i) = batch(@transient_fjrom_call, 0, {sel_method, i});
end


