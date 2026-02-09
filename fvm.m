function A = fvm(grid)

r = zeros(grid.numNodes,1); c = r; v = r;
qn = zeros(1,4);
entry = 1;

for E = 1:grid.numElements
    for sigma = 1:4
        switch sigma
            case 1
                neighbour = 2;
            case 2
                neighbour = 3;
            case 3
                neighbour = 4;
            case 4
                neighbour = 1;
        end
        %N = grid.N(subVolume,:,E);
        Nx = grid.Nx(sigma,:,E);
        Ny = grid.Ny(sigma,:,E);

        qn(1) = grid.D(E) * grid.L(E,sigma) * [Nx(1),Ny(1)] * grid.n(:,sigma,E);
        qn(2) = grid.D(E) * grid.L(E,sigma) * [Nx(2),Ny(2)] * grid.n(:,sigma,E);
        qn(3) = grid.D(E) * grid.L(E,sigma) * [Nx(3),Ny(3)] * grid.n(:,sigma,E);
        qn(4) = grid.D(E) * grid.L(E,sigma) * [Nx(4),Ny(4)] * grid.n(:,sigma,E);

        r(entry) = grid.elements(E,sigma); c(entry) = grid.elements(E,1); v(entry) = qn(1); entry = entry + 1;
        r(entry) = grid.elements(E,sigma); c(entry) = grid.elements(E,2); v(entry) = qn(2); entry = entry + 1;
        r(entry) = grid.elements(E,sigma); c(entry) = grid.elements(E,3); v(entry) = qn(3); entry = entry + 1;
        r(entry) = grid.elements(E,sigma); c(entry) = grid.elements(E,4); v(entry) = qn(4); entry = entry + 1;

        r(entry) = grid.elements(E,neighbour); c(entry) = grid.elements(E,1); v(entry) = -qn(1); entry = entry + 1;
        r(entry) = grid.elements(E,neighbour); c(entry) = grid.elements(E,2); v(entry) = -qn(2); entry = entry + 1;
        r(entry) = grid.elements(E,neighbour); c(entry) = grid.elements(E,3); v(entry) = -qn(3); entry = entry + 1;
        r(entry) = grid.elements(E,neighbour); c(entry) = grid.elements(E,4); v(entry) = -qn(4); entry = entry + 1;
    end
end

A = sparse(r,c,v,grid.numNodes,grid.numNodes);
A = A ./ grid.V;

end