clear all;

% Geometry of finite elements - triangulation, local global node numbering,
% coordinates, boundary nodes
[Coord,Elem,Nb,Db]=InitialMesh(1);

for nl=1:6
%% Edge-Node-Element Connections
[n2ed,ed2el]=edge(Elem,Coord);
%% Element Redrefine
[Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);

%%%%%%%%%%% find h %%%%%%%%%%%%%%
h(nl)=sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);

% No of degrees of freedom (initially solution at all the nodes are assumed as unknowns, dirichlet boundary conditions 
% if any are to be incorporated later on )
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes, unique(Db));

% Intializing the matrices
A=sparse(size(Coord,1),size(Coord,1)); % A is global stiffness matrix
b=sparse(size(Coord,1),1); % global load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly of A 
% stima is element stiffness matrices
for j=1:size(Elem,1),
    A(Elem(j,:),Elem(j,:))=A(Elem(j,:),Elem(j,:))+...
                                    stima(Coord(Elem(j,:),:));                                                                
end

% Assembly of load vector b
for j=1:size(Elem,1),
    b(Elem(j,:))=b(Elem(j,:))+det([1 1 1; Coord(Elem(j,:),:)'])*f(sum(Coord(Elem(j,:),:))/3)/6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann Conditions
if (~isempty(Nb))
   for j=1:size(Nb,1)
       b(Nb(j,:))= b(Nb(j,:))+norm(Coord(Nb(j,1),:)-Coord(Nb(j,2),:))*...
                                     u_N(Coord(Nb(j,1),:),Coord(Nb(j,2),:))/2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh=zeros(length(FullNodes),1);
% Dirichlet Conditions
if (~isempty(Db))
    Dbnodes=unique(Db);
    for j=1:size(Dbnodes,1)
       uh(Dbnodes(j),1)=ue(Coord(Dbnodes(j),:));
    end
end
b=b-A*uh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solving the linear system
uh(FreeNodes)=A(FreeNodes,FreeNodes)\b(FreeNodes);

% Exact solution at the nodes
u=u_nodes(Coord);

[L2e(nl),H1e(nl)]=Err(Coord, Elem, uh, u);

end

for ml=1:(nl-1)
    ocl2(ml)=log(L2e(ml)/L2e(ml+1))/log(h(ml)/h(ml+1))
    och1(ml)=log(H1e(ml)/H1e(ml+1))/log(h(ml)/h(ml+1))
end

 