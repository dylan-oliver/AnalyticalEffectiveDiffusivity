function coef = yWCoef(k,alpha,xv,yv,x,y)
%YWCOEF(k,alpha,xv,yv,x,y) returns the coefficient of the kth series term
%from the y-part of the series solution (resulting from the top boundary)
%in block (i,j) (with corresponding domain [xv(1),xv(2);yv(1),yv(2)]),
%evaluated at the point(s) (x,y).
x0 = xv(1); xL = xv(2);
y0 = yv(1); yW = yv(2);
lam = k * pi / (xL - x0);
coef = 2 * sin(lam*(x - x0)) .* (exp(-lam*(yW - y)) + exp(-lam*(yW + y - 2*y0))) / (k * pi * alpha * (1 - exp(-2*lam*(yW - y0))));
end