% If generating new data (generate = true), comment out lines 1-8  in run.m
% before running. Otherwise, ensure relevant data files are present in the
% data folder. Data should be generated once (generate = true) to create
% the necessary files.

clearvars
close all
clc

path = pwd();

warning('off')
spparms('piv_tol',0.5)

generate = false;

problemId = 1:9;

numTermsMarch = [2*(2:51)-3,2*52-4];

numTerms{1} = unique([10:5:400,numTermsMarch]); problem_x_ticks{1} = [1,100,200,300,400];
numTerms{2} = unique([10:5:400,numTermsMarch]); problem_x_ticks{2} = [1,100,200,300,400];
numTerms{3} = unique([10:5:400,numTermsMarch]); problem_x_ticks{3} = [1,100,200,300,400];
numTerms{4} = unique([10:5:300,numTermsMarch]); problem_x_ticks{4} = [1,100,200,300];
numTerms{5} = unique([10:5:200,numTermsMarch]); problem_x_ticks{5} = [1,50,100,150,200];
numTerms{6} = unique([10:2:100,numTermsMarch]); problem_x_ticks{6} = [1,20,40,60,80,100];
numTerms{7} = unique([10:2:70,numTermsMarch(1:35)]); problem_x_ticks{7} = [1,10,20,30,40,50,60,70];
numTerms{8} = unique([10:2:70,numTermsMarch(1:35)]); problem_x_ticks{8} = [1,10,20,30,40,50,60,70];
numTerms{9} = unique([10:2:52,numTermsMarch(1:26)]); problem_x_ticks{9} = [1,10,17,24,31,38,45,52];

if generate

    % Problems 1-9
    for id = problemId

        D_eff_analytical = zeros(2,2,length(numTerms{id}));

        timers_analytical = zeros(1,length(numTerms{id}));

        for n = 1:length(numTerms{id})

            clearvars -except path id n problemId numTerms numTermsMarch problem_x_ticks d Omega D_eff_analytical timers_analytical

            problem = id;
            S = numTerms{id}(n);

            analytical_method

            grid = store{1}.grid;
            D = store{1}.D;
            g = store{1}.g;
            alphay0 = store{1}.alphay0;
            alphayW = store{1}.alphayW;

            D_eff_analytical(:,:,n) = D_eff;

            timers_analytical(n) = timer_analytical;

        end

        clearvars -except path id n problemId numTerms numTermsMarch problem_x_ticks d Omega D_eff_analytical timers_analytical diffusivity_pattern grid g alphay0 alphayW

        problem = id;
        S = numTerms{id}(end);

        numerical_method

        D_eff_numerical = D_eff_fvm;

        timers_numerical = timer_fvm;

        save(['data/approximation_data_problem_',num2str(id),'.mat'],'D_eff_numerical','D_eff_analytical','timers_numerical','timers_analytical')
        save(['data/diffusivity_data_problem_',num2str(id),'.mat'],'D','diffusivity_pattern')
        save(['data/solution_data_problem_',num2str(id),'.mat'],'grid','g','alphay0','alphayW')

        disp(['Problem ',num2str(id),' complete...'])

    end

    % Problem 5 with maximised (rectangular) block size
    clearvars -except path problemId numTerms numTermsMarch problem_x_ticks
    close all
    clc

    problem = 10;
    S = 200;

    analytical_method

    numerical_method

    D_eff_numerical = D_eff_fvm;
    D_eff_analytical = D_eff;

    timers_numerical = timer_fvm;
    timers_analytical = timer_analytical;

    save('data/approximation_data_problem_5_rectangles.mat','D_eff_numerical','D_eff_analytical','timers_numerical','timers_analytical')
    save('data/diffusivity_data_problem_5_rectangles.mat','D','diffusivity_pattern')
    save('data/solution_data_problem_5_rectangles.mat','grid','g','alphay0','alphayW')

end

differences = cell(length(problemId));
timers = differences;
c1 = viridis(2);
c2 = plasma(2);
cols = [c1(1,:);c2(2,:)];
a_idx = [4,4,2,1,1,1,1,1,1];
n_idx = [6,3,2,2,2,1,1,1,1];
for id = problemId

    figure
    box on
    load(['data/diffusivity_data_problem_',num2str(id),'.mat'],'diffusivity_pattern')
    patch(diffusivity_pattern,'facecolor','flat','edgecolor','none')
    xlim([0,1])
    ylim([0,1])
    xticks([])
    xticklabels({})
    yticks([])
    yticklabels({})
    set(gcf,'position',[250,10,1200,1200])
    %set(gcf,'color','none')
    exportgraphics(gca,[path,'/Figures/diffusivity_pattern_',num2str(id),'.eps'],'ContentType','vector','BackgroundColor','none')

    load(['data/approximation_data_problem_',num2str(id),'.mat'],'D_eff_numerical','D_eff_analytical','timers_analytical')
    load(['data/approximation_data_problem_',num2str(id),'_march.mat'],'D_eff_numerical_march','D_eff_semianalytical_march','timers_march')

    createDiffusivityTex([path,'/Data'],['case_',num2str(id),'_analytical'],D_eff_analytical(:,:,end))
    createDiffusivityTex([path,'/Data'],['case_',num2str(id),'_numerical'],D_eff_numerical)

    if id == 6
        timers_analytical(end) = [];
        D_eff_analytical(:,:,end) = [];
        numTerms{id}(end) = [];
    end

    differences = abs(D_eff_analytical) - abs(D_eff_numerical);

    tmp = zeros(1,size(differences,3));
    norms_analytical = tmp;
    entry_differences = repmat(tmp,4,1);
    for n = 1:length(numTerms{id})

        tmp = differences(:,:,n);
        norms_analytical(n) = norm(tmp,'fro');
        entry_differences(:,n) = abs(tmp(:));

    end


    differences_march = abs(D_eff_semianalytical_march) - abs(D_eff_numerical_march);

    tmp = zeros(1,size(differences_march,3));
    norms_semianalytical_march = tmp;
    entry_differences_march = repmat(tmp,4,1);
    for n = 1:size(differences_march,3)

        tmp = differences_march(:,:,n);
        norms_semianalytical_march(n) = norm(tmp,'fro');
        entry_differences_march(:,n) = abs(tmp(:));

    end


    figure
    hold on
    box on
    plot(numTerms{id},norms_analytical,'color',cols(1,:),'linewidth',12)
    plot(numTermsMarch,norms_semianalytical_march,'color',cols(2,:),'linewidth',12)
    xticks(problem_x_ticks{id})
    xtickangle(0)
    if ismember(id,[7,8,9])
        tmp = xticklabels;
        tmp{end} = [tmp{end},'\phantom{0}'];
        xticklabels(tmp)
    end
    yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1])
    if id == 6
        xlim([min(numTerms{id}),100])
    else
        xlim([min(numTerms{id}),max(numTerms{id})])
    end
    ylim([1e-6,1e-0])
    set(gca,'yscale','log')
    xlabel('No. of series terms','interpreter','latex','fontsize',40)
    ylabel('$|| D_{\mathrm{eff}}^{(\mathrm{N})} - D_{\mathrm{eff}}^{(\mathrm{A})} ||_{\mathrm{F}}$','interpreter','latex','fontsize',40)
    set(gca,'ticklabelinterpreter','latex','fontsize',48,'xgrid','on','ygrid','on','zgrid','on','gridcolor','k','yscale','log')
    set(gcf,'position',[250,10,1200,1200])
    %set(gcf,'color','none')
    exportgraphics(gca,[path,'/Figures/D_eff_differences_',num2str(id),'.eps'],'ContentType','vector','BackgroundColor','none')


    figure
    hold on
    box on
    plot(timers_analytical(a_idx(id):end),norms_analytical(a_idx(id):end),'o','color',cols(1,:),'markerfacecolor',cols(1,:),'markersize',14)
    plot(timers_march(n_idx(id):end),norms_semianalytical_march(n_idx(id):end),'o','color',cols(2,:),'markerfacecolor',cols(2,:),'markersize',14)
    xticks([1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3])
    yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1])
    xlim([1e-3,1e3])
    ylim([1e-6,1e-0])
    set(gca,'yscale','log','xscale','log')
    xlabel('Time to compute $D_{\mathrm{eff}}^{(A)}$ (s)','interpreter','latex','fontsize',40)
    ylabel('$|| D_{\mathrm{eff}}^{(\mathrm{N})} - D_{\mathrm{eff}}^{(\mathrm{A})} ||_{\mathrm{F}}$','interpreter','latex','fontsize',40)
    set(gca,'ticklabelinterpreter','latex','fontsize',48,'xgrid','on','ygrid','on','zgrid','on','gridcolor','k','yscale','log')
    set(gcf,'position',[250,10,1200,1200])
    %set(gcf,'color','none')
    exportgraphics(gca,[path,'/Figures/D_eff_timers_',num2str(id),'.eps'],'ContentType','vector','BackgroundColor','none')


    [~,idx] = min(norms_analytical);
    createTex([path,'/Data'],['case_',num2str(id),'_analytical_norm_min_terms'],norms_analytical(idx))
    createTex([path,'/Data'],['case_',num2str(id),'_analytical_timer_min_terms'],timers_analytical(idx))


    [~,idx] = min(norms_semianalytical_march);
    createTex([path,'/Data'],['case_',num2str(id),'_semianalytical_norm_min_terms'],norms_semianalytical_march(idx))
    createTex([path,'/Data'],['case_',num2str(id),'_semianalytical_timer_min_terms'],timers_march(idx))

end

load('data/diffusivity_data_problem_5.mat','diffusivity_pattern')
figure
box on
patch(diffusivity_pattern,'facecolor','flat','edgecolor','k')
xlim([0,1])
ylim([0,1])
xticks([])
xticklabels({})
yticks([])
yticklabels({})
set(gcf,'position',[250,10,1200,1200])
%set(gcf,'color','none')
exportgraphics(gca,[path,'/Figures/diffusivity_pattern_5_comparison.eps'],'ContentType','vector','BackgroundColor','none')

load('data/diffusivity_data_problem_5_rectangles.mat','diffusivity_pattern')
figure
box on
patch(diffusivity_pattern,'facecolor','flat','edgecolor','k')
xlim([0,1])
ylim([0,1])
xticks([])
xticklabels({})
yticks([])
yticklabels({})
set(gcf,'position',[250,10,1200,1200])
%set(gcf,'color','none')
exportgraphics(gca,[path,'/Figures/diffusivity_pattern_5_rects.eps'],'ContentType','vector','BackgroundColor','none')

plot_density = 63;

for id = problemId

    load(['data/diffusivity_data_problem_',num2str(id),'.mat'],'D');
    D = D{1};
    D(D == 0.1) = 0;
    D(D == 1) = 1;
    image = logical(D);

    load(['data/solution_data_problem_',num2str(id),'.mat'],'grid','g','alphay0','alphayW')
    surface = @(x,y) genSol(grid,g,alphay0,alphayW,x,y);
    x = linspace(0,1,plot_density);
    y = linspace(0,1,plot_density);
    z = zeros(plot_density,plot_density);
    c = z;
    for i = 1:plot_density
        for j = 1:plot_density
            z(i,j) = surface(x(j),y(i));
            c(i,j) = imageDiffusivity(image,[0,1;0,1],[0.1,1],x(j),y(i));
        end
    end
    [x,y] = meshgrid(x,y);

    [faces,vertices,colours] = surf2patch(x,y,z,c);

    figure
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
    %set(gcf,'color','none')
    exportgraphics(gca,[path,'/Figures/solution_pattern_',num2str(id),'.eps'],'ContentType','vector','BackgroundColor','none')

end

% for id = problemId
% 
%     clearvars -except path problemId plot_density id
% 
%     load(['data/solution_data_problem_',num2str(id),'.mat'],'grid','g','alphay0','alphayW')
%     surface = @(x,y) genSol(grid,g,alphay0,alphayW,x,y);
%     x = linspace(0,1,plot_density);
%     y = linspace(0,1,plot_density);
%     z = zeros(plot_density,plot_density);
%     for i = 1:plot_density
%         for j = 1:plot_density
%             z(i,j) = surface(x(j),y(i));
%         end
%     end
%     [x,y] = meshgrid(x,y);
%     figure
%     surf(x,y,z,'facecolor','none','edgecolor','flat')
%     colormap(viridis(1000))
%     xlim([0,1])
%     ylim([0,1])
%     zlim([0,1])
%     xticks(0:0.2:1)
%     yticks(0:0.2:1)
%     zticks(0:0.2:1)
%     xlabel('$x$','interpreter','latex','fontsize',40)
%     ylabel('$y$','interpreter','latex','fontsize',40)
%     zlabel('$u^{(x)}(x,y)$','interpreter','latex','fontsize',40)
%     set(gca,'ticklabelinterpreter','latex','fontsize',48,'xgrid','on','ygrid','on','zgrid','on','gridcolor','k')
%     view([-60,10])
%     set(gca,'xgrid','on','ygrid','on','zgrid','on')
%     set(gcf,'position',[250,10,1200,1200])
%     exportgraphics(gca,[path,'/Figures/solution_pattern_',num2str(id),'.eps'],'ContentType','vector','BackgroundColor','none')
% 
% end



function createTex(path,name,value)

fid = fopen([path,'/',name,'.tex'],'w');

fprintf(fid,'%.4f',value);

fclose(fid);

end

function createDiffusivityTex(path,name,mat)

fid = fopen([path,'/',name,'.tex'],'w');

fprintf(fid,'\\begin{pmatrix} ');

if abs(mat(1,1)) < 5e-5
    fprintf(fid,'\\phantom{-}*');
else
    if mat(1,1) > 0
        fprintf(fid,'\\phantom{-}');
    end
    fprintf(fid,'%.4f',mat(1,1));
end

fprintf(fid,' & ');

if abs(mat(1,2)) < 5e-5
    fprintf(fid,'\\phantom{-}*');
else
    if mat(1,2) > 0
        fprintf(fid,'\\phantom{-}');
    end
    fprintf(fid,'%.4f',mat(1,2));
end

fprintf(fid,' \\\\ ');

if abs(mat(2,1)) < 5e-5
    fprintf(fid,'\\phantom{-}*');
else
    if mat(2,1) > 0
        fprintf(fid,'\\phantom{-}');
    end
    fprintf(fid,'%.4f',mat(2,1));
end

fprintf(fid,' & ');

if abs(mat(2,2)) < 5e-5
    fprintf(fid,'\\phantom{-}*');
else
    if mat(2,2) > 0
        fprintf(fid,'\\phantom{-}');
    end
    fprintf(fid,'%.4f',mat(2,2));
end

fprintf(fid,' \\end{pmatrix}');

fclose(fid);

end
