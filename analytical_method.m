%% Clear workspace
clearvars
close all
clc

%% Choose problem and corresponding parameters
problem = 1;

S = 3;

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

    case 2 % Centre square
        tmp = [
            d(1),d(1),d(1),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(1),d(1),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = D{1};

    case 3 % L-shape
        tmp = [
            d(1),d(1),d(1),d(1)
            d(1),d(2),d(2),d(1)
            d(1),d(1),d(2),d(1)
            d(1),d(1),d(2),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

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

    case 7 % Slice from image of geometry pattern
        load('data/imageSlice.mat')

        slice = slice(1:20,1:20);

        tmp = zeros(size(slice));

        tmp(slice) = d(1);
        tmp(~slice) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

    case 8 % Slice from image of geometry pattern
        load('data/imageSlice.mat')

        slice = slice(42:61,64:83);

        tmp = zeros(size(slice));

        tmp(slice) = d(1);
        tmp(~slice) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

    case 9 % Slice from image of geometry pattern
        load('data/imageSlice.mat')

        slice = slice(55:84,132:161);

        tmp = zeros(size(slice));

        tmp(slice) = d(1);
        tmp(~slice) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

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

    case 11 % Five vertical layers
        tmp = zeros(5,5);
        tmp(:,[1,3,5]) = d(1);
        tmp(:,[2,4]) = d(2);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

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

    case 13 % Inverted problem 1
        tmp = [
            d(2),d(1),d(2)
            d(1),d(2),d(1)
            d(2),d(1),d(2)
            ];

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

    case 14 % Checkerboard variant
        id = 12;
        tmp = load('data/checkerboard_investigation.mat','D');
        tmp = tmp.D(:,:,id);
        tmp(tmp == 1) = d(2);
        tmp(tmp == 0) = d(1);

        D{1} = flipud(tmp);

        D{2} = transpose(D{1});

    case 15 % Checkerboard 4x4
        tmp = [
            d(1),d(2),d(1),d(2)
            d(2),d(1),d(2),d(1)
            d(1),d(2),d(1),d(2)
            d(2),d(1),d(2),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = D{1};

    case 16 % Checkerboard 5x5
        tmp = [
            d(1),d(2),d(1),d(2),d(1)
            d(2),d(1),d(2),d(1),d(2)
            d(1),d(2),d(1),d(2),d(1)
            d(2),d(1),d(2),d(1),d(2)
            d(1),d(2),d(1),d(2),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = D{1};

    case 17 % Checkerboard 2x2
        tmp = [
            d(1),d(2)
            d(2),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = D{1};

    case 18
        tmp = [
            d(1),d(2),d(1)
            d(1),d(2),d(1)
            d(1),d(2),d(1)
            ];

        D{1} = flipud(tmp);

        D{2} = D{1};

end

tic

domain(:,:,1) = Omega;
domain(:,:,2) = Omega([2,1],:);

store = cell(1,2);

D_eff = zeros(2,2);

timer_part_1 = toc;

timer_part_2 = zeros(1,2);

for variant = 1:2

    tic

    [grid,right_interface_A0,right_interface_An,right_interface_Bn,top_interface_A0,top_interface_An,top_interface_Bn,b,alphay0,alphayW] = setConstants(domain(:,:,variant),D{variant},Sx,Sy);

    A = genA(grid,right_interface_A0,right_interface_An,right_interface_Bn,top_interface_A0,top_interface_An,top_interface_Bn);

    g = A \ b;

    D_eff(:,variant) = effectiveDiffusivity(grid,D{variant},g,alphay0,alphayW);

    timer_part_2(variant) = toc;

    store{variant}.grid = grid;
    store{variant}.D = D{variant};
    store{variant}.g = g;
    store{variant}.alphay0 = alphay0;
    store{variant}.alphayW = alphayW;

end

tic

D_eff([1,2],2) = D_eff([2,1],2);

timer_part_3 = toc;

timer_analytical = timer_part_1 + sum(timer_part_2) + timer_part_3;

%% Generate figures

% Diffusivity pattern for selected problem
cols = [
    0.336462660285714,0.336462660285714,0.336462660285714
    0.582485066142857,0.582485066142857,0.582485066142857
    ];

xSegments = zeros(4,store{1}.grid.L*store{1}.grid.W);
ySegments = zeros(4,store{1}.grid.L*store{1}.grid.W);
cSegments = zeros(store{1}.grid.L*store{1}.grid.W,3);
segmentIndex = 1:(4*store{1}.grid.L*store{1}.grid.W);
for i = 1:store{1}.grid.W

    for j = 1:store{1}.grid.L

        xSegments(:,j + (i-1)*store{1}.grid.L) = [store{1}.grid.xv(j);store{1}.grid.xv(j+1);store{1}.grid.xv(j+1);store{1}.grid.xv(j)];
        ySegments(:,j + (i-1)*store{1}.grid.L) = [store{1}.grid.yv(i);store{1}.grid.yv(i);store{1}.grid.yv(i+1);store{1}.grid.yv(i+1)];

        if D{1}(i,j) < max(d)
            cSegments(j + (i-1)*store{1}.grid.L,:) = cols(1,:);
        else
            cSegments(j + (i-1)*store{1}.grid.L,:) = cols(2,:);
        end

    end

end
diffusivity_pattern.faces = reshape(segmentIndex,4,store{1}.grid.L*store{1}.grid.W)';
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

% Increase to generate smoother 3D plots for x-/y-variant solutions below
plot_density = 65;

% x-variant solution
figure
image = D{1};
image(image == 0.1) = 0;
image(image == 1) = 1;
image = logical(image);
surface = @(x,y) genSol(store{1}.grid,store{1}.g,store{1}.alphay0,store{1}.alphayW,x,y);
x = linspace(0,1,plot_density);
y = linspace(0,1,plot_density);
z = zeros(plot_density,plot_density);
c = z;
for i = 1:plot_density
    for j = 1:plot_density
        z(i,j) = surface(x(j),y(i));
        c(i,j) = imageDiffusivity(image,Omega,d,x(j),y(i));
    end
end
[x,y] = meshgrid(x,y);
[faces,vertices,colours] = surf2patch(x,y,z,c);
patch('faces',faces,'vertices',vertices,'facevertexcdata',colours,'facecolor','none','edgecolor','flat','linewidth',5)
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
set(gcf,'position',[250,10,1200,1200])

% y-variant solution
figure
surface = @(x,y) genSol(store{1}.grid,store{1}.g,store{1}.alphay0,store{1}.alphayW,x,y);
x = linspace(0,1,plot_density);
y = linspace(0,1,plot_density);
z = zeros(plot_density,plot_density);
c = z;
image = transpose(image);
for i = 1:plot_density
    for j = 1:plot_density
        z(i,j) = surface(x(j),y(i));
        c(i,j) = imageDiffusivity(image,Omega,d,x(j),y(i));
    end
end
[x,y] = meshgrid(x,y);
[faces,vertices,colours] = surf2patch(x,y,z,c);
patch('faces',faces,'vertices',vertices,'facevertexcdata',colours,'facecolor','none','edgecolor','flat','linewidth',5)
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
set(gcf,'position',[250,10,1200,1200])
