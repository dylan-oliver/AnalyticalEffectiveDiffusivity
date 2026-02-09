function grad = effectiveFVM(u,grid)
numNodes = grid.numNodes;
numElements = grid.numElements;
grad = zeros(numNodes,2);
for E = 1:numElements
    element = grid.elements(E,:);
    vals = u(element);
    D = grid.D(E);
    for subVol = 1:4
        Nx = grid.Nx_centroid(subVol,:,E);
        Ny = grid.Ny_centroid(subVol,:,E);
        grad(element(subVol),1) = grad(element(subVol),1) + grid.a(subVol,E) * D * Nx * vals;
        grad(element(subVol),2) = grad(element(subVol),2) + grid.a(subVol,E) * D * Ny * vals;
    end
end
grad = sum(grad,1);
grad = grad(:);
end
