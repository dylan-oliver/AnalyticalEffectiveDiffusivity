% %% Clear workspace
% clearvars
% close all
% clc
% 
% %% Choose problem and corresponding parameters
% problem = 3;
% 
% S = 200;

Sx = S;
Sy = S;

d(1) = 1;
d(2) = 0.1;

Omega = [
    0,1
    0,1
    ];

switch problem
    case 1 % Checkerboard
        tmp = [
            d(1),d(2),d(1),d(2)
            d(2),d(1),d(2),d(1)
            d(1),d(2),d(1),d(2)
            d(2),d(1),d(2),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = D{1};

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 2 % Centre square
        tmp = [
            d(1),d(1),d(1),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(1),d(1),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = D{1};

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 3 % L-shape
        tmp = [
            d(1),d(1),d(1),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(1),d(2),d(1)
            d(1),d(1),d(2),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 4 % Cross
        tmp = [
            d(1),d(1),d(2),d(1),d(1)
            d(2),d(2),d(2),d(2),d(2)
            d(1),d(1),d(2),d(1),d(1)
            d(1),d(1),d(2),d(1),d(1)
            d(1),d(1),d(2),d(1),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 5 % Centre/corner cross pattern
        tmp = [
            d(2),d(2),d(1),d(1),d(1),d(1),d(2),d(2)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(2),d(2),d(2),d(2),d(1),d(1)
            d(1),d(1),d(2),d(2),d(2),d(2),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(2),d(2),d(1),d(1),d(1),d(1),d(2),d(2)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 6 % Pattern from N. March (2021)
        tmp1 = [
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(2),d(2),d(2),d(2),d(2),d(1)
            d(1),d(2),d(2),d(2),d(2),d(2),d(2),d(1)
            d(1),d(2),d(2),d(2),d(2),d(2),d(2),d(1)
            d(1),d(1),d(2),d(2),d(2),d(2),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(2),d(2),d(2),d(1)
            ];

        tmp2 = [
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(2),d(2),d(2),d(2),d(2),d(1),d(1)
            d(2),d(2),d(2),d(2),d(2),d(2),d(1),d(1)
            d(1),d(2),d(2),d(2),d(2),d(2),d(1),d(1)
            d(1),d(1),d(2),d(2),d(2),d(2),d(1),d(1)
            d(1),d(1),d(1),d(2),d(2),d(1),d(1),d(1)
            d(1),d(1),d(2),d(2),d(2),d(1),d(1),d(1)
            ];

        tmp3 = [
            d(1),d(1),d(2),d(2),d(2),d(2),d(2),d(2)
            d(1),d(1),d(2),d(2),d(2),d(2),d(2),d(2)
            d(1),d(1),d(2),d(2),d(2),d(2),d(2),d(2)
            d(1),d(1),d(2),d(2),d(2),d(2),d(2),d(2)
            d(1),d(1),d(1),d(2),d(2),d(2),d(2),d(2)
            d(1),d(1),d(1),d(1),d(2),d(2),d(2),d(2)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(2)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            ];

        tmp4 = [
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(2),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(2),d(2),d(2),d(2),d(2),d(1),d(1),d(1)
            d(2),d(2),d(2),d(2),d(2),d(2),d(1),d(1)
            d(2),d(2),d(2),d(2),d(2),d(2),d(1),d(1)
            d(2),d(2),d(2),d(2),d(2),d(1),d(1),d(1)
            d(2),d(2),d(2),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            ];

        tmp = [tmp1,tmp2;tmp3,tmp4];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 7 % Slice from image of geometry pattern
        load('data/imageSlice.mat')

        slice = slice(1:20,1:20);

        tmp = zeros(size(slice));

        tmp(slice) = d(1);
        tmp(~slice) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = flipud(slice);
        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 8 % Slice from image of geometry pattern
        load('data/imageSlice.mat')

        slice = slice(42:61,64:83);

        tmp = zeros(size(slice));

        tmp(slice) = d(1);
        tmp(~slice) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = flipud(slice);
        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 9 % Slice from image of geometry pattern
        load('data/imageSlice.mat')

        slice = slice(55:84,132:161);

        tmp = zeros(size(slice));

        tmp(slice) = d(1);
        tmp(~slice) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = flipud(slice);
        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 10 % Centre/corner cross pattern, maximised block size
        tmp = [
            d(2),d(1),d(1),d(2)
            d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1)
            d(2),d(1),d(1),d(2)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 11 % Five vertical layers
        tmp = zeros(5,5);
        tmp(:,[1,3,5]) = d(1);
        tmp(:,[2,4]) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 12 % Three horizontal layers
        tmp = [
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(2),d(2),d(2),d(2),d(2),d(2),d(2),d(2)
            d(2),d(2),d(2),d(2),d(2),d(2),d(2),d(2)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            d(1),d(1),d(1),d(1),d(1),d(1),d(1),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 13 % Inverted problem 1
        tmp = [
            d(2),d(1),d(2)
            d(1),d(2),d(1)
            d(2),d(1),d(2)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 14 % Checkerboard variant
        id = 12;
        tmp = load('data/checkerboard_investigation.mat','D');
        tmp = tmp.D(:,:,id);
        image = flipud(logical(tmp));
        tmp(tmp == 1) = d(2);
        tmp(tmp == 0) = d(1);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 15 % Checkerboard 4x4
        tmp = [
            d(1),d(2),d(1),d(2)
            d(2),d(1),d(2),d(1)
            d(1),d(2),d(1),d(2)
            d(2),d(1),d(2),d(1)
            ];

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D{1} = flipud(tmp);

        D{2} = D{1};

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

    case 16 % Checkerboard 5x5
        tmp = [
            d(1),d(2),d(1),d(2),d(1)
            d(2),d(1),d(2),d(1),d(2)
            d(1),d(2),d(1),d(2),d(1)
            d(2),d(1),d(2),d(1),d(2)
            d(1),d(2),d(1),d(2),d(1)
            ];

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D{1} = flipud(tmp);

        D{2} = D{1};

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);


    case 17 % Checkerboard 2x2
        tmp = [
            d(1),d(2)
            d(2),d(1)
            ];

        image = tmp;
        image(image < 1) = 0;
        image(image > 0) = 1;
        image = flipud(logical(image));

        D{1} = flipud(tmp);

        D{2} = D{1};

        D_fvm = @(x,y) imageDiffusivity(image,Omega,d,x,y);

end

tic

L = size(D{1},2);
W = size(D{1},1);

x0 = Omega(1,1); xL = Omega(1,2); xNumNodes = Sy * L + 1;
y0 = Omega(2,1); yH = Omega(2,2); yNumNodes = Sx * W + 1;

num_fvm_nodes = xNumNodes * yNumNodes;

conditions(1,:) = [1,1,2,2];
conditions(2,:) = [2,2,1,1];

g_fvm{1}{1} = @(x,y) x0 + zeros(size(x));
g_fvm{1}{2} = @(x,y) xL + zeros(size(x));
g_fvm{1}{3} = @(x,y) zeros(size(y));
g_fvm{1}{4} = @(x,y) zeros(size(y));

g_fvm{2}{1} = @(x,y) zeros(size(x));
g_fvm{2}{2} = @(x,y) zeros(size(x));
g_fvm{2}{3} = @(x,y) y0 + zeros(size(y));
g_fvm{2}{4} = @(x,y) yH + zeros(size(y));

grid_fvm = genGrid(D_fvm,x0,xL,y0,yH,xNumNodes,yNumNodes);
A_fvm_unmodified = fvm(grid_fvm);

sol_fvm = zeros(num_fvm_nodes,2);
D_eff_fvm = zeros(2,2);

for variant = 1:2

    [A_fvm,b_fvm] = setBoundary(grid_fvm,conditions(variant,:),g_fvm{variant},A_fvm_unmodified);

    sol_fvm(:,variant) = A_fvm \ b_fvm;

    D_eff_fvm(:,variant) = effectiveFVM(sol_fvm(:,variant),grid_fvm);

end

timer_fvm = toc;

%% Generate figures

% Diffusivity pattern for selected problem
figure
cols = [
    0.336462660285714,0.336462660285714,0.336462660285714
    0.582485066142857,0.582485066142857,0.582485066142857
    ];

xv = Omega(1,1):(Omega(1,2) - Omega(1,1)) / L:Omega(1,2);
yv = Omega(2,1):(Omega(2,2) - Omega(2,1)) / W:Omega(2,2);

xSegments = zeros(4,L*W);
ySegments = zeros(4,L*W);
cSegments = zeros(L*W,3);
segmentIndex = 1:(4*L*W);
for i = 1:W

    for j = 1:L

        xSegments(:,j + (i-1)*L) = [xv(j);xv(j+1);xv(j+1);xv(j)];
        ySegments(:,j + (i-1)*L) = [yv(i);yv(i);yv(i+1);yv(i+1)];

        if D{1}(i,j) < max(d)
            cSegments(j + (i-1)*L,:) = cols(1,:);
        else
            cSegments(j + (i-1)*L,:) = cols(2,:);
        end

    end

end
diffusivity_pattern.faces = reshape(segmentIndex,4,L*W)';
diffusivity_pattern.vertices = [xSegments(:),ySegments(:)];
diffusivity_pattern.facevertexcdata = cSegments;
patch(diffusivity_pattern,'facecolor','flat','edgecolor','none')
xlim([0,1])
ylim([0,1])
xticks([])
xticklabels({})
yticks([])
yticklabels({})
set(gcf,'position',[250,100,1200,800])

% x-variant solution
figure
c = zeros(length(grid_fvm.nodes(:,1)),1);
for i = 1:length(grid_fvm.nodes(:,1))
    x = grid_fvm.nodes(i,1);
    y = grid_fvm.nodes(i,2);
    c(i) = D_fvm(x,y);
end
patch('faces',grid_fvm.elements,'vertices',[grid_fvm.nodes,sol_fvm(:,1)],'facevertexcdata',c,'facecolor','none','edgecolor','flat')
colormap(viridis(2))
xlim([0,1])
ylim([0,1])
zlim([0,1])
xticks(0:0.2:1)
yticks(0:0.2:1)
zticks(0:0.2:1)
xlabel('$x$','interpreter','latex','fontsize',40)
ylabel('$y$','interpreter','latex','fontsize',40)
zlabel('$u^{(x)}(x,y)$','interpreter','latex','fontsize',40)
set(gca,'ticklabelinterpreter','latex','fontsize',48,'xgrid','on','ygrid','on','zgrid','on','gridcolor','k')
view([-60,10])
set(gca,'xgrid','on','ygrid','on','zgrid','on')
set(gcf,'position',[250,100,1200,800])

% y-variant solution
figure
c = zeros(length(grid_fvm.nodes(:,1)),1);
for i = 1:length(grid_fvm.nodes(:,1))
    x = grid_fvm.nodes(i,1);
    y = grid_fvm.nodes(i,2);
    c(i) = D_fvm(x,y);
end
patch('faces',grid_fvm.elements,'vertices',[grid_fvm.nodes,sol_fvm(:,2)],'facevertexcdata',c,'facecolor','none','edgecolor','flat')
colormap(viridis(2))
xlim([0,1])
ylim([0,1])
zlim([0,1])
xticks(0:0.2:1)
yticks(0:0.2:1)
zticks(0:0.2:1)
xlabel('$x$','interpreter','latex','fontsize',40)
ylabel('$y$','interpreter','latex','fontsize',40)
zlabel('$u^{(x)}(x,y)$','interpreter','latex','fontsize',40)
set(gca,'ticklabelinterpreter','latex','fontsize',48,'xgrid','on','ygrid','on','zgrid','on','gridcolor','k')
view([-60,10])
set(gca,'xgrid','on','ygrid','on','zgrid','on')
set(gcf,'position',[250,100,1200,800])
