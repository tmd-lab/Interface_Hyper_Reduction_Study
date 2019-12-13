function [OutputNodes,OutputQuad,OutputTri,NodeCorrelation] = CreatOutputFiles(filename,NewNodes,NewElements,OldNodes,varargin)
% CreatOutputFiles: Creates output files of reduced meshes
% INPUTS:
%   filename : Filename of output
%   NewNodes : Array with nodes of reduced mesh
%   NewElements : New element with reduced nodes
%   OldNodes : Array with nodes of original mesh
%   varargin    : Node IDS for reduced mesh (defalt 1:No. Nodes)
% OUTPUTS: File containing
%   OutputNodes: Nodes [idx,x,y,z]
%   OutputQuad: List with quadrilateral elements
%   OutputTri: List with triangular elements
%   NodeCorrelation: list with idex of original nodes located in reduced
%   mesh's elements

filename=[filename,num2str(length(NewElements)),'Elements_'];
filext={'Nodes.dat','QuadEle.dat','TriEle.dat','OldNewRelation.dat'};
if isempty(varargin)
    IDnewnodes=[1:length(NewNodes)]';
else
    IDnewnodes=cell2mat(varargin);
    if ~all(size(IDnewnodes)==[length(NewNodes),1])
        error('Size of ID Vector is not correct');
    end
end    


OutputNodes=[IDnewnodes,NewNodes(:,:)];  % Nodes Format: ID, x , y , z
% csvwrite([filename,filext{1}],OutputNodes)

% Initializing
OutputQuadTemp=zeros(length(NewElements),5);
OutputTriTemp=zeros(length(NewElements),4);
OutputQuadRefTemp=cell(length(NewElements),1);
OutputTriRefTemp=cell(length(NewElements),1);
nquad=0;
ntri=0;
sortidx4=[1:4 1:4];
sortidx3=[1:3 1:3];
for ii=1:length(NewElements)
    if length(NewElements{ii,1})==4
        nquad=nquad+1;
        % Sort Nodes: (Node 1 to Node 2 greatest x change)
        [~,iidx]=max(diff([OutputNodes(NewElements{ii,1},2);OutputNodes(NewElements{ii,1}(1),2)]));
        OutputQuadTemp(nquad,2:end)=OutputNodes(NewElements{ii,1}(sortidx4(iidx:iidx+3)),1);
        OutputQuadRefTemp{nquad,1}=NewElements{ii,2};
    else
        ntri=ntri+1;
        % Sort Nodes: (Node 1 to Node 2 greatest x change)
        [~,iidx]=max(diff([OutputNodes(NewElements{ii,1},2);OutputNodes(NewElements{ii,1}(1),2)]));     
        OutputTriTemp(ntri,2:end)=OutputNodes(NewElements{ii,1}(sortidx3(iidx:iidx+2)),1);
        OutputTriRefTemp{ntri,1}=NewElements{ii,2};
    end
end
% Add Element ID
IDQuads=[1:nquad]';
IDTri=[nquad+1:nquad+ntri]';
OutputQuad=[IDQuads,OutputQuadTemp(1:nquad,2:end)];
OutputTri=[IDTri,OutputTriTemp(1:ntri,2:end)];

% Output Element List - 
% Format: Element ID [1:No.Elements], NodeID1, Node ID 2 , Node ID3, ..]
% csvwrite([filename,filext{2}],OutputQuad) 
% csvwrite([filename,filext{3}],OutputTri)

% Old-New Node correlation
Temp=cell(size(OldNodes,1),1);
for k=1:nquad
    for i=1:length(OutputQuadRefTemp{k,1})
        Temp{OutputQuadRefTemp{k,1}(i),1}=[Temp{OutputQuadRefTemp{k,1}(i),1},IDQuads(k)];
    end
end
for k=1:ntri
    for i=1:length(OutputTriRefTemp{k,1})
        Temp{OutputTriRefTemp{k,1}(i),1}=[Temp{OutputTriRefTemp{k,1}(i),1},IDTri(k)];
    end
end
%Determine Number of Colums for Element ID Output 
L=cellfun(@length,Temp);
if min(L)==0
    warning('Output is wrong');
end
Output4EleID=zeros(size(OldNodes,1),max(L));
for k=1:length(Output4EleID)
    Output4EleID(k,1:L(k))=Temp{k,1};
end
NodeCorrelation=[OldNodes,Output4EleID];
% Format ID , x , y , z , Elem ID1 , Elem ID2 , ...
% (First for Colums are from original Nodeset)
% csvwrite([filename,filext{4}],NodeCorrelation) 
filename=filename(1:end-1);

% % Plot each patch seperatly to check if correct nodes were found.
% figure
% for kk=1:nquad
% clf
% plot(OldNodes(:,2),OldNodes(:,3),'r*');
% hold on
% polyElem=polyshape(OutputNodes(OutputQuad(kk,2:end),2),OutputNodes(OutputQuad(kk,2:end),3));
% plot(polyElem);
% plot(OldNodes(OutputQuadRefTemp{kk},2),OldNodes(OutputQuadRefTemp{kk},3),'bo');
% axis equal
% pause(1);
% end
% 
% for kk=1:ntri
% clf
% plot(OldNodes(:,2),OldNodes(:,3),'r*');
% hold on
% polyElem=polyshape(OutputNodes(OutputTri(kk,2:end),2),OutputNodes(OutputTri(kk,2:end),3));
% plot(polyElem);
% plot(OldNodes(OutputTriRefTemp{kk},2),OldNodes(OutputTriRefTemp{kk},3),'bo');
% axis equal
% pause(1);
% end

save(filename,'OutputNodes','OutputQuad','OutputTri','NodeCorrelation');
disp(['The Output for ',num2str(nquad),' quadrilateral and ',num2str(ntri),' triangular Elements was generated!']);
end

