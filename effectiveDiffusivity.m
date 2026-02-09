function A = effectiveDiffusivity(grid,D,g,alphay0,alphayH)


L = grid.L;
W = grid.W;
Nx = grid.Nx;
Ny = grid.Ny;
index_right = grid.index_right;
index_top = grid.index_top;
xv = grid.xv;
yv = grid.yv;
A = zeros(2,1);

for i = 1:W

    for j = 1:L

        if j == 1
            coefA0 = xv(j) * (yv(i+1) - yv(i));
        else
            coefA0 = g(index_right(i,j-1));
        end
        A(1) = A(1) - D(i,j) * coefA0;

        if j == L
            coefB0 = xv(j+1) * (yv(i+1) - yv(i));
        else
            coefB0 = g(index_right(i,j));
        end
        A(1) = A(1) + D(i,j) * coefB0;

        for k = 1:Nx

            if j == 1
                coefAn = 0;
            else
                coefAn = g(index_right(i,j-1) + k);
            end
            A(2) = A(2) + 2 * D(i,j) * (1 + (-1)^(k+1)) * (2 * exp(-k*pi*(xv(j+1) - xv(j))/(yv(i+1) - yv(i))) - exp(-2*k*pi*(xv(j+1) - xv(j))/(yv(i+1) - yv(i))) - 1) / (k * pi * (1 - exp(-2*k*pi*(xv(j+1) - xv(j))/(yv(i+1) - yv(i))))) * coefAn;

            if j == L
                coefBn = 0;
            else
                coefBn = g(index_right(i,j) + k);
            end
            A(2) = A(2) + 2 * D(i,j) * (1 + (-1)^(k+1)) * (2 * exp(-k*pi*(xv(j+1) - xv(j))/(yv(i+1) - yv(i))) - exp(-2*k*pi*(xv(j+1) - xv(j))/(yv(i+1) - yv(i))) - 1) / (k * pi * (1 - exp(-2*k*pi*(xv(j+1) - xv(j))/(yv(i+1) - yv(i))))) * coefBn;

        end

        for k = 1:Ny

            if i == 1
                coefCn = 0;
            else
                coefCn = g(index_top(i-1,j) + k);
            end
            A(2) = A(2) - 2 * D(i,j) * (xv(j+1) - xv(j)) * (1 + (-1)^(k+1)) * (2 * exp(-k*pi*(yv(i+1) - yv(i))/(xv(j+1) - xv(j))) - exp(-2*k*pi*(yv(i+1) - yv(i))/(xv(j+1) - xv(j))) - 1) / (k^2 * pi^2 * alphay0(i,j) * (1 - exp(-2*k*pi*(yv(i+1) - yv(i))/(xv(j+1) - xv(j))))) * coefCn;

            if i == W
                coefDn = 0;
            else
                coefDn = g(index_top(i,j) + k);
            end
            A(2) = A(2) - 2 * D(i,j) * (xv(j+1) - xv(j)) * (1 + (-1)^(k+1)) * (2 * exp(-k*pi*(yv(i+1) - yv(i))/(xv(j+1) - xv(j))) - exp(-2*k*pi*(yv(i+1) - yv(i))/(xv(j+1) - xv(j))) - 1) / (k^2 * pi^2 * alphayH(i,j) * (1 - exp(-2*k*pi*(yv(i+1) - yv(i))/(xv(j+1) - xv(j))))) * coefDn;

        end

    end

end

end