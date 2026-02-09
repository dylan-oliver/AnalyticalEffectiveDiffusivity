% If generating new data (generate = true), comment out lines 1-8  in run.m
% before running. Otherwise, ensure relevant data files are present in the
% data folder. Data should be generated once (generate = true) to create
% the necessary files.

clearvars
close all
clc

path = pwd();

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

yMax = [30,30,30,40,100,100,100,100,200];

differences = cell(length(problemId));
timers = differences;
c1 = viridis(2);
c2 = plasma(2);
cols = [c1(1,:);c2(2,:)];
a_idx = [4,4,2,1,1,1,1,1,1];
n_idx = [6,3,2,2,2,1,1,1,1];
for id = problemId

    load(['data/approximation_data_problem_',num2str(id),'.mat'],'D_eff_numerical','D_eff_analytical','timers_analytical')
    load(['data/approximation_data_problem_',num2str(id),'_march.mat'],'D_eff_numerical_march','D_eff_semianalytical_march','timers_march')

    if id == 6
        timers_analytical(end) = [];
        D_eff_analytical(:,:,end) = [];
        numTerms{id}(end) = [];
    end

    figure
    hold on
    box on
    plot(numTerms{id}(a_idx(id):end),timers_analytical(a_idx(id):end),'o','color',cols(1,:),'markerfacecolor',cols(1,:),'markersize',14)
    plot(numTermsMarch(n_idx(id):end),timers_march(n_idx(id):end),'o','color',cols(2,:),'markerfacecolor',cols(2,:),'markersize',14)
    xticks(problem_x_ticks{id})
    yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1])
    if id == 6
        xlim([min(numTerms{id}),100])
    else
        xlim([min(numTerms{id}),max(numTerms{id})])
    end
    ylim([0,yMax(id)])
    set(gca,'yscale','log')
    xlabel('No. of series terms','interpreter','latex','fontsize',40)
    ylabel('Time to compute $D_{\mathrm{eff}}^{(A)}$ (s)','interpreter','latex','fontsize',40)
    set(gca,'ticklabelinterpreter','latex','fontsize',48,'xgrid','on','ygrid','on','zgrid','on','gridcolor','k','yscale','log')
    set(gcf,'position',[250,10,1200,1200])

end
