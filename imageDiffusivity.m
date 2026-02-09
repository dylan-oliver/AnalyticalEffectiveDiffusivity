function D = imageDiffusivity(image,domain,d,x,y)

L = size(image,2);
W = size(image,1);

x_index = L*x / (domain(1,2) - domain(1,1)); x_index = floor(x_index) + 1; x_index = min(L,x_index);
y_index = W*y / (domain(2,2) - domain(2,1)); y_index = floor(y_index) + 1; y_index = min(W,y_index);

if image(y_index,x_index)
    D = d(1);
else
    D = d(2);
end

end