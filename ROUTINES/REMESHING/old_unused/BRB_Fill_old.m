
load('Input_Greg')
d=Elements_idx;
b=inp;
[unique_elements,unique_nodes,corner_nodes]=BRB_Fill(Element,eList)
function [unique_elements,unique_nodes,corner_nodes]= BRB_Fill(Element,eList)
unique_elements=unique(eList);%remove duplicates in input vector
node_set=[];
for c=1:length(unique_elements)
    
    node_set=vertcat(node_set,Element{unique_elements(c),1}(:)); %creating node set from all input elements
end
unique_nodes=unique(node_set);   
e=1;
while e<=length(Element)
    if length(intersect(Element{e,1},unique_nodes))==3%add element if there are 3 shared nodes /w node set and test element
        unique_elements=horzcat(unique_elements,e);
        node_set=vertcat(node_set,Element{e,1}');
        unique_nodes=unique(node_set);
        e=0;
    end
    e=e+1;
end
corner_nodes=[];
for x=1:length(node_set)
    y=find(node_set(x)==node_set);
    if length(y)==1
        corner_nodes=vertcat(corner_nodes,node_set(x));
    end
end
end


