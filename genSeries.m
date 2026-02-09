function series_solution = genSeries(grid,g,alphay0,alphayW,x,y)

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

left_linear = zeros(2,1);
right_linear = left_linear;
left = zeros(2,Nx);
right = left;
bottom = zeros(2,Ny);
top = bottom;

% Linear term from left-hand boundary
if j == 1
    coefA0 = xv(j) * (yv(i+1) - yv(i));
else
    coefA0 = g(index_right(i,j-1));
end
left_linear(1) = coefA0;
left_linear(2) = -(xTmp - xv(j+1)) / ((xv(j+1) - xv(j)) * (yv(i+1) - yv(i)));

% Linear term from right-hand boundary
if j == L
    coefB0 = xv(j+1) * (yv(i+1) - yv(i));
else
    coefB0 = g(index_right(i,j));
end
right_linear(1) = coefB0;
right_linear(2) = (xTmp - xv(j)) / ((xv(j+1) - xv(j)) * (yv(i+1) - yv(i)));


for k = 1:Nx
    % Series solution from left-hand boundary
    if j == 1
        coefAn = 0;
    else
        coefAn = g(index_right(i,j-1) + k);
    end
    left(1,k) = coefAn;
    left(2,k) = x0Coef(k,xv(j:j+1),yv(i:i+1),xTmp,yTmp);

    % Series solution from right-hand boundary
    if j == L
        coefBn = 0;
    else
        coefBn = g(index_right(i,j) + k);
    end
    right(1,k) = coefBn;
    right(2,k) = xLCoef(k,xv(j:j+1),yv(i:i+1),xTmp,yTmp);
end

for k = 1:Ny
    % Series solution from bottom boundary
    if i == 1
        coefCn = 0;
    else
        coefCn = g(index_top(i-1,j) + k);
    end
    bottom(1,k) = coefCn;
    bottom(2,k) = y0Coef(k,alphay0(i,j),xv(j:j+1),yv(i:i+1),xTmp,yTmp);

    % Series solution from top boundary
    if i == W
        coefDn = 0;
    else
        coefDn = g(index_top(i,j) + k);
    end
    top(1,k) = coefDn;
    top(2,k) = yWCoef(k,alphayW(i,j),xv(j:j+1),yv(i:i+1),xTmp,yTmp);
end

series_solution.left_linear = left_linear;
series_solution.right_linear = right_linear;
series_solution.left = left;
series_solution.right = right;
series_solution.bottom = bottom;
series_solution.top = top;

end