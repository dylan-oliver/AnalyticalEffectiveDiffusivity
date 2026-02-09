function solTmp = genSol(grid,g,alphay0,alphayW,x,y)

L = grid.L;
W = grid.W;
Nx = grid.Nx;
Ny = grid.Ny;

xv = grid.xv;
yv = grid.yv;

index_right = grid.index_right;
index_top = grid.index_top;

i = W*y / (yv(end) - yv(1)); i = floor(i) + 1; i = min(W,i);
j = L*x / (xv(end) - xv(1)); j = floor(j) + 1; j = min(L,j);

xTmp = x;
yTmp = y;
solTmp = 0*x;

% Linear term from left-hand boundary
if j == 1
    coefA0 = xv(j) * (yv(i+1) - yv(i));
else
    coefA0 = g(index_right(i,j-1));
end
solTmp = solTmp - (xTmp - xv(j+1)) / ((xv(j+1) - xv(j)) * (yv(i+1) - yv(i))) * coefA0;

% Linear term from right-hand boundary
if j == L
    coefB0 = xv(j+1) * (yv(i+1) - yv(i));
else
    coefB0 = g(index_right(i,j));
end
solTmp = solTmp + (xTmp - xv(j)) / ((xv(j+1) - xv(j)) * (yv(i+1) - yv(i))) * coefB0;


for k = 1:Nx
    % Series solution from left-hand boundary
    if j == 1
        coefAn = 0;
    else
        coefAn = g(index_right(i,j-1) + k);
    end
    solTmp = solTmp + x0Coef(k,xv(j:j+1),yv(i:i+1),xTmp,yTmp) * coefAn;

    % Series solution from right-hand boundary
    if j == L
        coefBn = 0;
    else
        coefBn = g(index_right(i,j) + k);
    end
    solTmp = solTmp + xLCoef(k,xv(j:j+1),yv(i:i+1),xTmp,yTmp) * coefBn;
end

for k = 1:Ny
    % Series solution from bottom boundary
    if i == 1
        coefCn = 0;
    else
        coefCn = g(index_top(i-1,j) + k);
    end
    solTmp = solTmp + y0Coef(k,alphay0(i,j),xv(j:j+1),yv(i:i+1),xTmp,yTmp) * coefCn;

    % Series solution from top boundary
    if i == W
        coefDn = 0;
    else
        coefDn = g(index_top(i,j) + k);
    end
    solTmp = solTmp + yWCoef(k,alphayW(i,j),xv(j:j+1),yv(i:i+1),xTmp,yTmp) * coefDn;
end

end