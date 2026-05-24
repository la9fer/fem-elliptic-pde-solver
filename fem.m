%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM for  - div . (u_x,u_y)   = f    in  \Omega             %
%                            u = 0   on  \partial\Omega_D   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
format long;

% Geometry of finite elements - triangulation, local global node numbering,
% coordinates, boundary nodes

[Coord,Elem,Nb,Db]=InitialMesh(1);
%figure
%triplot(Elem,Coord(:,1),Coord(:,2))

maxlevel=5;
L2e=zeros(1,maxlevel); H1e=zeros(1,maxlevel);
for i=1:maxlevel
  
%% Edge-Node-Element Connections
[n2ed,ed2el]=edge(Elem,Coord); %nodes2edges, edge2element
%% Element Redrefine
[Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);
%  figure
%  triplot(Elem,Coord(:,1),Coord(:,2))

di(i)=1/(2^i);

% No of degrees of freedom (initially solution at all the nodes
%are assumed as unknowns, dirichlet boundary conditions 
% if any are to be incorporated later on )
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes, unique(Db));

% Intializing the matrices
A=sparse(size(Coord,1),size(Coord,1)); % A is global stiffness matrix
b=sparse(size(Coord,1),1); % global load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly of A 
% stima is element stiffness matrices
for j=1:size(Elem,1)
    A(Elem(j,:),Elem(j,:))=A(Elem(j,:),Elem(j,:))+...
                                    stima(Coord(Elem(j,:),:));                                                                
end

% Assembly of load vector b
for j=1:size(Elem,1)
    b(Elem(j,:))=b(Elem(j,:))+det([1 1 1; Coord(Elem(j,:),:)'])*f(sum(Coord(Elem(j,:),:))/3)/6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solving the linear system
uh(FreeNodes)=A(FreeNodes,FreeNodes)\b(FreeNodes);
% Exact solution at the nodes
u=u_nodes(Coord);
[L2e(i), H1e(i)]=Err(Coord, Elem, uh, u); %error computation

end
 
 ocuh1=zeros(1,maxlevel-1); ocul2=zeros(1,maxlevel-1);
 for jj=1:(maxlevel-1) %convergence rate
    ocuh1(jj)=log(H1e(jj)/H1e(jj+1))/log(di(jj)/di(jj+1));
    ocul2(jj)=log(L2e(jj)/L2e(jj+1))/log(di(jj)/di(jj+1));
  end
 
 ocuh1
 ocul2

% % Display the computed solution
figure
show(Coord,Elem,full(uh),full(u))




