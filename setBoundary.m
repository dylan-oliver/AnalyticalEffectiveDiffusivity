function [A,b] = setBoundary(grid,conditions,g,A)

b = zeros(grid.numNodes,1);
dirichlet = false(grid.numNodes,1);
for E = 1:grid.numEdges
    edge = grid.boundary(E,:);
    if conditions(edge(4)) == 1
        b(edge(1:2)) = g{edge(4)}(grid.nodes(edge(1:2),1),grid.nodes(edge(1:2),2));
        dirichlet(edge(1:2)) = true;
    end
end

A(dirichlet,:) = 0;
I = speye(size(A));
I(~dirichlet,:) = 0;
A = I + A;

end