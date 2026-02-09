function grid = genGrid(D,x0,xL,y0,yW,xNumNodes,yNumNodes)
x = linspace(x0,xL,xNumNodes);

y = linspace(y0,yW,yNumNodes);

numNodes = xNumNodes * yNumNodes;
numElements = (xNumNodes - 1) * (yNumNodes - 1);
numEdges = 2*(xNumNodes + yNumNodes) - 4;

nodes = zeros(numNodes,2);
elements = zeros(numElements,4);
boundary = zeros(numEdges,4);

DElement = zeros(numElements,1);

index = @(i,j) i + (j - 1)*yNumNodes;
elementId = 1;
edgeId = 1;
cols = 1:4;
for j = 1:xNumNodes-1
    for i = 1:yNumNodes-1
        idx = [index(i,j), index(i,j+1), index(i+1,j+1), index(i+1,j)];
        nodes(idx,:) = [x(j), y(i);x(j+1), y(i);x(j+1), y(i+1);x(j), y(i+1)];
        elements(elementId,cols) = idx;

        DElement(elementId) = D((x(j) + x(j+1))/2, (y(i) + y(i+1))/2);

        if any(nodes(idx,1) == x0)
            boundary(edgeId,:) = [idx(4), idx(1), elementId, 1];
            edgeId = edgeId + 1;
        end

        if any(nodes(idx,1) == xL)
            boundary(edgeId,:) = [idx(2), idx(3), elementId, 2];
            edgeId = edgeId + 1;
        end

        if any(nodes(idx,2) == y0)
            boundary(edgeId,:) = [idx(1), idx(2), elementId, 3];
            edgeId = edgeId + 1;
        end

        if any(nodes(idx,2) == yW)
            boundary(edgeId,:) = [idx(3), idx(4), elementId, 4];
            edgeId = edgeId + 1;
        end

        elementId = elementId + 1;
    end
end

grid.numNodes = numNodes;
grid.numElements = numElements;
grid.numEdges = numEdges;
grid.nodes = nodes;
grid.elements = elements;
grid.boundary = sortrows(boundary,4);
grid.D = DElement;

grid = setVals(grid);

end
