addpath C:\Users\Dave\Downloads\matgraph % adds matgraph to the matlab script.
maxdegree = input('enter the integer maximum degree of the tree T \n');
%maxdegree = 3;
n = input('enter an integer,n, the number of nodes in T \n');
%n=20;
%alpha = input('enter the degree of preferential atachment 0 < \alpha < 2 \n');
alpha = 1;

% eventually we will add options for inputing other parameters such as
% tortuosity and radii.
%  
g = graph;
add(g,1,2);
add(g,2,3);

%Now we make the tree g.

 %Initially we set g to be a graph with three  nodes 1,2 and 3 such that {1,2} and {2,3} are 
 %connected by a single edge.
%what follows is the set up with preferential attachment.

j=1;
while j <= n-3,
    v(1) = (deg(g,1))^(alpha);  %v is a vector.  First entry v(1) is the degree of node 1 taken to the power of apha
    
    for k = 2:j+2,              %this runs from 2 to n-1 eventually
        v(k) = (deg(g,k))^alpha + v(k-1);
    end
    
    max = v(j+2);
    
    r = rand;   %chooses a number uniformly between 0 and 1
    y = r*max;
    
    w = v;
    w(j+3) = y;
    x = sort(w);

l=1;
    while x(l) < y,
    l = l+1;
    end
    if l ~= 1;
        if deg(g,l) < maxdegree,
            add(g,l,j+3) 
            j=j+1;
        end
    else
        if deg(g,l) < (maxdegree - 1),
            add(g,l,j+3) 
            j=j+1;
        end    
    end    
end

nauty(g,'roxdraft.dre');
%at this point it might be useful to save the graph in nauty format, open
%dreadnaut, calculate the automorphsim group and see what happens.

clf
ndraw(g); %draws the tree g.  We should design an option for this to be switched on/off earlier on.

A = double(matrix(g));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            %a rough approximation of cts time
%%%%%%%%%Now we concentrate on the concentration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = zeros(1,100);

S = 95;
for t = 1:100;             %a rough approximation of cts time
    if ne(g) ~= 0,
        for i=1:n-1,
            for j=i+1:n,
                if A(i,j)==1,
                    s=randi(100);
                    if s>99,
                        delete(g,i,j); %llok up matgraph delete edge i j
                    end
                end
            end
        end
        if isconnected(g)==0,  %if graph g is not connected
            B = components(g);  %produces a family of sets , each
                   %set a connected component of 
            for m=3:n,
                con = size(B(1));
                conner = size(B(m));
                if conner(2) ~= con(2),  %If the component contaning node m  is not the same size as the component contating node 1 then we delete edges incident to m.
                    secdegree = deg(g,m);
                    neigh=neighbors(g,m);
                    for q=1:secdegree,  
                        delete(g,m,neigh(q));
                    end
                else   %If the component contaning the root and the component containg node m are the same size.  Check if they are the same component.  If not delete edges incident to m.
                    if B(m) ~= B(1),         
                        thirddegree = deg(g,m);
                        neighs=neighbors(g,m);
                        for s=1:thirddegree,
                            delete(g,m,neighs(s));
                        end
                    end
                end    
            end
        end
    F(t) = n - 1 - ne(g);%This measures the total number of edges that have been removed at each time t
    A=double(matrix(g));
    end
end
