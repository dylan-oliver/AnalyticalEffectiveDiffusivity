function coef = xLCoef(k,xv,yv,x,y)
%X0COEF(k,xv,yv,x,y) returns the coefficient of the kth series term from
%the x-part of the series solution (resulting from the right boundary) in
%block (i,j) (with corresponding domain [xv(1),xv(2);yv(1),yv(2)]),
%evaluated at the point(s) (x,y).
x0 = xv(1); xL = xv(2);
y0 = yv(1); yW = yv(2);
lam = k * pi / (yW - y0);
coef = 2 * (exp(-lam*(xL - x)) - exp(-lam*(xL + x - 2*x0))) .* cos(lam*(y - y0)) / ((yW - y0) * (1 - exp(-2*lam*(xL - x0))));
end