function grid = setVals(grid)
nodes = grid.nodes;
elements = grid.elements;

numNodes = grid.numNodes;
numElements = grid.numElements;

L = zeros(numElements,4);
n = zeros(2,4,numElements);
a = zeros(4,numElements);
V = zeros(numNodes,1);
C = zeros(4,4,numElements);
N = zeros(4,4,numElements);
Nx = zeros(4,4,numElements);
Ny = zeros(4,4,numElements);
Nx_centroid = zeros(4,4,numElements);
Ny_centroid = zeros(4,4,numElements);

mid = zeros(4,2);
sv_centroid = zeros(4,2);
onesPad = [1;1;1;1];
zeroPad = [0;0;0;0];
for E = 1:numElements
    % Extract grid element and corresponding vertices
    square = elements(E,:);
    vertex = nodes(square(1:4),:);

    % Element mid-point
    centre = mean(vertex,1);

    % Define element side mid-points
    s1 = (vertex(2,:) + vertex(1,:))/2;
    s2 = (vertex(3,:) + vertex(2,:))/2;
    s3 = (vertex(4,:) + vertex(3,:))/2;
    s4 = (vertex(1,:) + vertex(4,:))/2;

    % Define sub-volume edges
    e1 = centre - s1;
    e2 = centre - s2;
    e3 = centre - s3;
    e4 = centre - s4;

    % Determine sub-volume edge mid-points
    mid(1,:) = s1 + e1/2;
    mid(2,:) = s2 + e2/2;
    mid(3,:) = s3 + e3/2;
    mid(4,:) = s4 + e4/2;

    % Determine sub-volume mid-points
    sv_centroid(1,:) = [mid(4,1),mid(1,2)];
    sv_centroid(2,:) = [mid(2,1),mid(1,2)];
    sv_centroid(3,:) = [mid(2,1),mid(3,2)];
    sv_centroid(4,:) = [mid(4,1),mid(3,2)];

    % Determine element edge magnitudes
    L(E,1) = norm(e1);
    L(E,2) = norm(e2);
    L(E,3) = norm(e3);
    L(E,4) = norm(e4);

    % Define and normalise corresponding edge normals
    n(:,1,E) = [e1(2);-e1(1)]/L(E,1);
    n(:,2,E) = [e2(2);-e2(1)]/L(E,2);
    n(:,3,E) = [e3(2);-e3(1)]/L(E,3);
    n(:,4,E) = [e4(2);-e4(1)]/L(E,4);

    % Determine sub-volume areas
    a(1,E) = L(E,1) * L(E,4);
    a(2,E) = L(E,1) * L(E,2);
    a(3,E) = L(E,2) * L(E,3);
    a(4,E) = L(E,3) * L(E,4);

    % Contribute sub-volume areas to overall control volume area
    V(square(1:4)) = V(square(1:4)) + a(:,E);

    % Evaluating interpolant s(x,y) = a*x + b*y + c*x*y + d in each
    % element, where x = {x1,x2,x3,x4} and y = {y1,y2,y3,y4} are the (x,y)
    % points stored in 'vertex' and {u1,u2,u3,u4} are the corresponding
    % solution values at those vertices, yields the linear system:
    % [ x1, y1, x1*y1, 1 ] [a]   [u1]
    % [ x2, y2, x2*y2, 1 ] [b]   [u2]
    % [ x3, y3, x3*y3, 1 ] [c] = [u3]
    % [ x4, y4, x4*y4, 1 ] [d]   [u4]

    % C is the transposed cofactor matrix of the above coefficient matrix,
    % and detC is its determinant.
    C(:,:,E) = [vertex(:,1),vertex(:,2),vertex(:,1).*vertex(:,2),onesPad];

    % After solving for the interpolant coefficients a, b, c, and d (see
    % above) in terms of solutions u1, u2, u3, and u4 at corresponding
    % element vertices, the resulting function is rearranged, forming the
    % shape function s(u1,u2,u3,u4) = N1*u1 + N2*u2 + N3*u3 + N4*u4.
    % The rows of N, Nx, and Ny correspond to the shape function evaluated
    % at the mid-point of each sub-volume edge.
    N(:,:,E) = [mid(:,1), mid(:,2), mid(:,1).*mid(:,2), onesPad] / C(:,:,E);
    Nx(:,:,E) = [onesPad, zeroPad, mid(:,2), zeroPad] / C(:,:,E);
    Ny(:,:,E) = [zeroPad, onesPad, mid(:,1), zeroPad] / C(:,:,E);

    % Also evaluate shape function at sub-volume mid-points for
    % approximation of effective diffusivity
    Nx_centroid(:,:,E) = [onesPad,zeroPad,sv_centroid(:,2),zeroPad] / C(:,:,E);
    Ny_centroid(:,:,E) = [zeroPad,onesPad,sv_centroid(:,1),zeroPad] / C(:,:,E);
end

grid.L = L;
grid.n = n;
grid.a = a;
grid.V = V;
grid.N = N;
grid.Nx = Nx;
grid.Ny = Ny;
grid.Nx_centroid = Nx_centroid;
grid.Ny_centroid = Ny_centroid;
grid.C = C;

end
