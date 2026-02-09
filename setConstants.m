function [grid,right_interface_A0,right_interface_An,right_interface_Bn,top_interface_A0,top_interface_An,top_interface_Bn,b,alphay0,alphayW] = setConstants(domain,D,Nx,Ny,xv,yv)
%SETCONSTANTS(domain,D,Nx,Ny,xv,yv) defines necessary values utilised by
%analytical solution (series solution coefficients, boundary coefficients,
%etc.)
%
% Input:
% domain - (2,2) array defined: [x0,xH;y0,yV], geometric properties of a
%          given domain.
% D      - (V,H) array, defines a pixelated heterogeneous domain (for
%          example, D=[1,0,1;0,1,0;1,0,1] defines a 'checkerboard' pattern
%          domain where material 1 is assigned to values where D=1, and
%          material 2 is assigned to values where D=0).
% Nx     - Number of terms/coefficients to define for the x-series part of
%          the solution. 
% Ny     - Number of terms/coefficients to define for the y-series part of
%          the solution. 
% xv     - (optional) (1,:) array of values corresponding to the x-location
%          of block/pixel interfaces/boundaries.
% yv     - (optional) (1,:) array of values corresponding to the y-location
%          of block/pixel interfaces/boundaries.
%
% Output:
% grid               - Structure containing information regarding the
%                      values associated with the geometric domain (domain
%                      lengths, interface locations, number of unknowns
%                      resulting from series solution, etc. See end of
%                      file.)
% right_interface_A0 - (V,H-1,4) array. Contains constant terms from x-part
%                      of series solution corresponding to each homogeneous
%                      block.
% right_interface_An - (Nx,Nx+1,V,H-1,4) array. Contains series
%                      coefficients terms from x-part of series solution
%                      corresponding to each homogeneous block.  
% right_interface_Bn - (Ny,Nx+1,V,H-1,4) array. Contains series
%                      coefficients terms from y-part of series solution 
%                      corresponding to each homogeneous block.
% top_interface_A0   - (Ny,V-1,H,4) array. Contains constant terms from
%                      x-part of series solution corresponding to each
%                      homogeneous block.
% top_interface_An   - (Nx,Ny,V-1,H,4) array. Contains series coefficients
%                      terms from x-part of series solution corresponding
%                      to each homogeneous block.
% top_interface_Bn   - (Ny,Ny,V-1,H,4) array. Contains series coefficients
%                      terms from y-part of series solution corresponding
%                      to each homogeneous block.
% b                  - (V*(H-1)*(Nx+1) + H*(V-1)*Ny,1) array. Contains
%                      constant terms from application of interface
%                      conditions to the series solutions of each block.
% alphay0            - (V,H) array. Contains either 1 or Dij corresponding
%                      to the condition enforced at the bottom side of
%                      block (i,j).
% alphayW            - (V,H) array. Contains either 1 or Dij corresponding
%                      to the condition enforced at the top side of block
%                      (i,j).

H = size(D,2); % Domain horizontal length
V = size(D,1); % Domain vertical length

%Indices right: from ((j - 1 + (i-1)*(L-1))*(N+1) + 1) to ((j + (i-1)*(L-1))*(N+1))
index_right = @(i,j) (j - 1 + (i-1)*(H-1))*(Nx+1) + 1;

%Indices top: from ((i - 1 + (j-1)*(W-1))*N + W*(L-1)*(N+1) + 1) to ((i + (j-1)*(W-1))*N + W*(L-1)*(N+1))
index_top = @(i,j) (i - 1 + (j-1)*(V-1))*Ny + V*(H-1)*(Nx+1);

total_unknowns = V*(H-1)*(Nx+1) + H*(V-1)*Ny;

if nargin < 5
    xv = domain(1,1):(domain(1,2) - domain(1,1)) / H:domain(1,2);
    yv = domain(2,1):(domain(2,2) - domain(2,1)) / V:domain(2,2);
end

% Right-hand side vector after gathering like-term coefficients when 
% applying interface conditions to series solutions in each block (all
% known constant values are pushed to RHS and are contained in b).
b = zeros(total_unknowns,1);

% Alpha is defined according to the location of a given 'pixel' or block 
% (see equations 30 and 31 in article)
alphay0 = zeros(V,H);
alphayW = zeros(V,H);
for i = 1:V

    for j = 1:H

        if i == 1
            alphay0(i,j) = 1;
        else
            alphay0(i,j) = D(i,j);
        end

        if i == V
            alphayW(i,j) = 1;
        else
            alphayW(i,j) = D(i,j);
        end

    end

end

% Set constants for x-variant solution
row_index_right = @(i,j) (j - 1 + (i-1)*(H-1))*(Nx+1);

% Series coefficients for x-variant solution
right_interface_A0 = zeros(V,H-1,4);
right_interface_An = zeros(Nx,Nx+1,V,H-1,4);
right_interface_Bn = zeros(Ny,Nx+1,V,H-1,4);

for i = 1:V

    for j = 1:H-1

        eta = (1:Nx) * pi / (yv(i+1) - yv(i));
        mu = (1:Ny) * pi / (xv(j+1) - xv(j));

        coef0 = 1 / ((xv(j+1) - xv(j)) * (yv(i+1) - yv(i)));
        coefx = 2 ./ ((yv(i+1) - yv(i)) * (1 - exp(-2*eta*(xv(j+1) - xv(j)))));
        coefy = 2 ./ ((xv(j+1) - xv(j)) * mu .* (1 - exp(-2*mu*(yv(i+1) - yv(i)))));

        x = xv(j+1);
        y = linspace(yv(i),yv(i+1),Nx+1);

        for k = 1:Nx

            right_interface_An(k,:,i,j,1) = -D(i,j) * coefx(k) * eta(k) * (exp(-eta(k)*(2*xv(j+1) - x - xv(j))) + exp(-eta(k)*(x - xv(j)))) .* cos(eta(k)*(y - yv(i)));
            right_interface_An(k,:,i,j,2) = D(i,j) * coefx(k) * eta(k) * (exp(-eta(k)*(xv(j+1) - x)) + exp(-eta(k)*(xv(j+1) + x - 2*xv(j)))) .* cos(eta(k)*(y - yv(i)));

            right_interface_An(k,:,i,j,3) = D(i,j+1) * coefx(k) * eta(k) * (exp(-eta(k)*(2*xv(j+2) - x - xv(j+1))) + exp(-eta(k)*(x - xv(j+1)))) .* cos(eta(k)*(y - yv(i)));
            right_interface_An(k,:,i,j,4) = -D(i,j+1) * coefx(k) * eta(k) * (exp(-eta(k)*(xv(j+2) - x)) + exp(-eta(k)*(xv(j+2) + x - 2*xv(j+1)))) .* cos(eta(k)*(y - yv(i)));

        end

        for k = 1:Ny

            right_interface_Bn(k,:,i,j,1) = -D(i,j) * coefy(k) * mu(k) * cos(mu(k)*(x - xv(j))) .* (exp(-mu(k)*(2*yv(i+1) - y - yv(i))) + exp(-mu(k)*(y - yv(i)))) / alphay0(i,j);
            right_interface_Bn(k,:,i,j,2) = D(i,j) * coefy(k) * mu(k) * cos(mu(k)*(x - xv(j))) .* (exp(-mu(k)*(yv(i+1) - y)) + exp(-mu(k)*(yv(i+1) + y - 2*yv(i)))) / alphayW(i,j);

            right_interface_Bn(k,:,i,j,3) = D(i,j+1) * coefy(k) * mu(k) * cos(mu(k)*(x - xv(j+1))) .* (exp(-mu(k)*(2*yv(i+1) - y - yv(i))) + exp(-mu(k)*(y - yv(i)))) / alphay0(i,j+1);
            right_interface_Bn(k,:,i,j,4) = -D(i,j+1) * coefy(k) * mu(k) * cos(mu(k)*(x - xv(j+1))) .* (exp(-mu(k)*(yv(i+1) - y)) + exp(-mu(k)*(yv(i+1) + y - 2*yv(i)))) / alphayW(i,j+1);

        end

        idx = row_index_right(i,j) + (1:Nx+1);

        % A0(i,j) values
        if j == 1
            b(idx) = b(idx) + D(i,j) * xv(j)*(yv(i+1) - yv(i)) * coef0;
        else
            right_interface_A0(i,j,1) = -D(i,j) * coef0;
        end

        % B0(i,j) values
        right_interface_A0(i,j,2) = D(i,j) * coef0;


        % A0(i,j+1) values
        right_interface_A0(i,j,3) = D(i,j+1) * coef0;

        % B0(i,j+1) values
        if j == H-1
            b(idx) = b(idx) + D(i,j+1) * xv(j+2)*(yv(i+1) - yv(i)) * coef0;
        else
            right_interface_A0(i,j,4) = -D(i,j+1) * coef0;
        end

    end

end

% Set constants for y-variant solution
row_index_top = @(i,j) (i - 1 + (j-1)*(V-1))*Ny + (H-1)*V*(Nx+1);

% Series coefficients for y-variant solution
top_interface_A0 = zeros(Ny,V-1,H,4);
top_interface_An = zeros(Nx,Ny,V-1,H,4);
top_interface_Bn = zeros(Ny,Ny,V-1,H,4);

for j = 1:H

    for i = 1:V-1

        eta = (1:Nx) * pi / (yv(i+1) - yv(i));
        mu = (1:Ny) * pi / (xv(j+1) - xv(j));

        coef0 = 1 / ((xv(j+1) - xv(j)) * (yv(i+1) - yv(i)));
        coefx = 2 ./ ((yv(i+1) - yv(i)) * (1 - exp(-2*eta*(xv(j+1) - xv(j)))));
        coefy = 2 ./ ((xv(j+1) - xv(j)) * mu .* (1 - exp(-2*mu*(yv(i+1) - yv(i)))));

        % Offset node on right-hand boundary (j==L) to avoid singular 
        % coefficient matrix
        if j == H
            x = linspace(xv(j),xv(j+1),Ny+1); x(end) = x(end) - 1e-14; x(1) = [];
        else
            x = linspace(xv(j),xv(j+1),Ny+1); x(1) = [];
        end
        y = yv(i+1);

        for k = 1:Nx

            top_interface_An(k,:,i,j,1) = -coefx(k) * (exp(-eta(k)*(2*xv(j+1) - x - xv(j))) - exp(-eta(k)*(x - xv(j)))) .* cos(eta(k)*(y - yv(i)));
            top_interface_An(k,:,i,j,2) = coefx(k) * (exp(-eta(k)*(xv(j+1) - x)) - exp(-eta(k)*(xv(j+1) + x - 2*xv(j)))) .* cos(eta(k)*(y - yv(i)));

            top_interface_An(k,:,i,j,3) = coefx(k) * (exp(-eta(k)*(2*xv(j+1) - x - xv(j))) - exp(-eta(k)*(x - xv(j)))) .* cos(eta(k)*(y - yv(i+1)));
            top_interface_An(k,:,i,j,4) = -coefx(k) * (exp(-eta(k)*(xv(j+1) - x)) - exp(-eta(k)*(xv(j+1) + x - 2*xv(j)))) .* cos(eta(k)*(y - yv(i+1)));

        end

        for k = 1:Ny

            top_interface_Bn(k,:,i,j,1) = -coefy(k) * sin(mu(k)*(x - xv(j))) .* (exp(-mu(k)*(2*yv(i+1) - y - yv(i))) + exp(-mu(k)*(y - yv(i)))) / alphay0(i,j);
            top_interface_Bn(k,:,i,j,2) = coefy(k) * sin(mu(k)*(x - xv(j))) .* (exp(-mu(k)*(yv(i+1) - y)) + exp(-mu(k)*(yv(i+1) + y - 2*yv(i)))) / alphayW(i,j);

            top_interface_Bn(k,:,i,j,3) = coefy(k) * sin(mu(k)*(x - xv(j))) .* (exp(-mu(k)*(2*yv(i+2) - y - yv(i+1))) + exp(-mu(k)*(y - yv(i+1)))) / alphay0(i+1,j);
            top_interface_Bn(k,:,i,j,4) = -coefy(k) * sin(mu(k)*(x - xv(j))) .* (exp(-mu(k)*(yv(i+2) - y)) + exp(-mu(k)*(yv(i+2) + y - 2*yv(i+1)))) / alphayW(i+1,j);

        end

        idx = row_index_top(i,j) + (1:Ny);

        % A0(i,j) values
        if j == 1
            b(idx) = b(idx) + xv(j) * (yv(i+1) - yv(i)) * (x(:) - xv(j+1)) * coef0;
        else
            top_interface_A0(:,i,j,1) = -(x(:) - xv(j+1)) * coef0;
        end

        % B0(i,j) values
        if j == H
            b(idx) = b(idx) - xv(j+1) * (yv(i+1) - yv(i)) * (x(:) - xv(j)) * coef0;
        else
            top_interface_A0(:,i,j,2) = (x(:) - xv(j)) * coef0;
        end

        % A0(i+1,j) values
        if j == 1
            b(idx) = b(idx) - xv(j) * (yv(i+2) - yv(i+1)) * (x(:) - xv(j+1)) * coef0;
        else
            top_interface_A0(:,i,j,3) = (x(:) - xv(j+1)) * coef0;
        end

        % B0(i+1,j) values
        if j == H
            b(idx) = b(idx) + xv(j+1) * (yv(i+2) - yv(i+1)) * (x(:) - xv(j)) * coef0;
        else
            top_interface_A0(:,i,j,4) = -(x(:) - xv(j)) * coef0;
        end

    end

end

grid.L = H;
grid.W = V;

grid.Nx = Nx;
grid.Ny = Ny;

grid.index_right = index_right;
grid.index_top = index_top;

grid.xv = xv;
grid.yv = yv;

grid.total_unknowns = total_unknowns;

end
