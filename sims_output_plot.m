% Most up to date script for plotting simulation outputs

blue    = '#0072BD';
orange  = '#D95319';
purple  = '#7E2F8E';
yellow  = '#EDB120';
green   = '#77AC30';
cyan    = '#4DBEEE';
colres  = 64;
col     = [max(min(3/2-abs(3-(1:colres)*4/colres),1),0)' max(min(3/2-abs(2-(1:colres)*4/colres),1),0)' max(min(3/2-abs(1-(1:colres)*4/colres),1),0)'];
iternum = 5e2;

%% ADULTS
ITERATIONNUMBER = 1; % CHANGE ACCORDING TO NUMBER OF ITERATIONS
pname = input('enter directory path name: ','s')
zj = -1;
za = -1;

for SAVINGPLOTS=1:ITERATIONNUMBER
    try
        fname = ['iter_' num2str(SAVINGPLOTS) '_data.mat'];
        iter = importdata(fullfile(pname,fname));

        time = 1:iternum;
        juvenilephenos = iter(2:iternum+1,2);
        adultphenos = iter(2:iternum+1,4);
        juveniles = iter(2:iternum+1,1);
        adults = iter(2:iternum+1,3);

        for h = 1:iternum
            if any(juveniles{h}==0)
                ext  = find(juveniles{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end

            if any(adults{h}==0)
                ext  = find(adults{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end
        end

        % Time series of lineages scaled by population size
        figure(2)
        clf()
        set(gcf,'Position',[100 100 700 700])
        iter = iternum;
        for i=1:iternum
            x_axis = NaN(iternum,1);
            y_axis = NaN(iternum,1);
            x_axis(1:numel(juvenilephenos{i})) = juvenilephenos{i};
            y_axis(1:numel(adultphenos{i})) = adultphenos{i};
            z_axis = zeros(iternum,1);
            z_axis(1:numel(z_axis)) = i;
            sizes  = nan(iternum,1);
            sizes(1:numel(adults{i})) = sqrt(adults{i}/sum(adults{i}))*2e2;
            sizes = sqrt(sizes)*3e1;
            scatter3(x_axis,y_axis,z_axis,sizes,'MarkerFaceColor',col(ceil(i/iter*colres),:),'MarkerEdgeColor','k')
            hold on
        end
        %ax=gca;
        %ax.FontSize = 16;
        axis([-0.01 1 -0.01 1 0 iternum])
        xticks([0:0.5:1])
        yticks([0:0.5:1])
        zticks([0:100:iternum])
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',24,'linewidth',7)
        xlabel('\Theta_J')
        ylabel('\Theta_A')
        zlabel('Time')
        %legend('Lineages over time')
        sgtitle([' '])%Adult dominated. Time series of lineages scaled by population size. zj = ' num2str(zj) ', za = ' num2str(za)])
        view(0,90)
        pbaspect([1 1 1])
        figname = ['iter_' num2str(SAVINGPLOTS) '_lineages.fig'];
        saveas(gcf,fullfile(pname,figname))
        figname2 = ['iter_' num2str(SAVINGPLOTS) '_lineages.pdf'];
        print(fullfile(pname,figname2),'-bestfit','-dpdf')
    catch
    end
end

%% PRODUCE SIMULATION FIGURES FOR ANIMATION
SAVINGPLOTS = 5
fname = ['iter_' num2str(SAVINGPLOTS) '_data.mat'];
        iter = importdata(fullfile(pname,fname));

        time = 1:iternum;
        juvenilephenos = iter(2:iternum+1,2);
        adultphenos = iter(2:iternum+1,4);
        juveniles = iter(2:iternum+1,1);
        adults = iter(2:iternum+1,3);

        for h = 1:iternum
            if any(juveniles{h}==0)
                ext  = find(juveniles{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end

            if any(adults{h}==0)
                ext  = find(adults{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end
        end
figure(2)
iter = iternum;
j=1;
for i=1:iternum
    clf()
    x_axis = NaN(iternum,1);
    y_axis = NaN(iternum,1);
    x_axis(1:numel(juvenilephenos{i})) = juvenilephenos{i};
    y_axis(1:numel(adultphenos{i})) = adultphenos{i};
    z_axis = zeros(iternum,1);
    z_axis(1:numel(z_axis)) = i;
    sizes  = nan(iternum,1);
    sizes(1:numel(adults{i})) = sqrt(adults{i}/sum(adults{i}))*2e2;
    sizes = sqrt(sizes)*3e1;
    scatter3(x_axis,y_axis,z_axis,sizes,'MarkerFaceColor',col(ceil(i/iter*colres),:),'MarkerEdgeColor','k')
    set(gcf,'Position',[100 100 700 700])
    axis([-0.01 1 -0.01 1 0 iternum])
    xticks([0:0.5:1])
    yticks([0:0.5:1])
    zticks([0:100:iternum])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',24,'linewidth',7)
    xlabel('\Theta_J')
    ylabel('\Theta_A')
    zlabel('Time')
    %legend('Lineages over time')
    sgtitle([' '])%Adult dominated. Time series of lineages scaled by population size. zj = ' num2str(zj) ', za = ' num2str(za)])
    view(0,90)
    pbaspect([1 1 1])
    %figname = ['iter_' num2str(SAVINGPLOTS) '_lineages_' num2str(j) '.pdf'];
    %print(fullfile(pname,figname),'-bestfit','-dpdf')
    figname = ['iter_' num2str(SAVINGPLOTS) '_lineages_' num2str(j) '.png'];
    saveas(gcf,fullfile(pname,figname))
    j=j+1;
end

%% PRODUCE GIF FROM SIMULATION DATA

SAVINGPLOTS = 5
fname = ['iter_' num2str(SAVINGPLOTS) '_data.mat'];
        iter = importdata(fullfile(pname,fname));

        time = 1:iternum;
        juvenilephenos = iter(2:iternum+1,2);
        adultphenos = iter(2:iternum+1,4);
        juveniles = iter(2:iternum+1,1);
        adults = iter(2:iternum+1,3);

        for h = 1:iternum
            if any(juveniles{h}==0)
                ext  = find(juveniles{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end

            if any(adults{h}==0)
                ext  = find(adults{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end
        end
obj = figure(1)
removeToolbarExplorationButtons(obj)
iter = iternum;
j=1;
gifFile = 'trimorphism.gif';
exportgraphics(obj, gifFile);
for i=1:iternum
    clf()
    x_axis = NaN(iternum,1);
    y_axis = NaN(iternum,1);
    x_axis(1:numel(juvenilephenos{i})) = juvenilephenos{i};
    y_axis(1:numel(adultphenos{i})) = adultphenos{i};
    z_axis = zeros(iternum,1);
    z_axis(1:numel(z_axis)) = i;
    sizes  = nan(iternum,1);
    sizes(1:numel(adults{i})) = sqrt(adults{i}/sum(adults{i}))*2e2;
    sizes = sqrt(sizes)*3e1;
    scatter3(x_axis,y_axis,z_axis,sizes,'MarkerFaceColor',col(ceil(i/iter*colres),:),'MarkerEdgeColor','k')
    %set(gcf,'Position',[100 100 700 700])
    axis([-0.01 1 -0.01 1 0 iternum])
    xticks([0:0.5:1])
    yticks([0:0.5:1])
    zticks([0:100:iternum])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',24,'linewidth',7)
    xlabel('\Theta_J')
    ylabel('\Theta_A')
    zlabel('Time')
    %legend('Lineages over time')
    sgtitle([' '])%Adult dominated. Time series of lineages scaled by population size. zj = ' num2str(zj) ', za = ' num2str(za)])
    view(0,90)
    pbaspect([1 1 1])
    %figname = ['iter_' num2str(SAVINGPLOTS) '_lineages_' num2str(j) '.pdf'];
    %print(fullfile(pname,figname),'-bestfit','-dpdf')
    if ~mod(j,3)
        exportgraphics(obj, gifFile, Append=true);   
    else
    end
    j=j+1;
end



%% JUVENILES
ITERATIONNUMBER = 10; % CHANGE ACCORDING TO NUMBER OF ITERATIONS
pname = input('enter directory path name: ','s')
% zj = -0.8;
% za = -0.2;

for SAVINGPLOTS=1:ITERATIONNUMBER
    try
        fname = ['iter_' num2str(SAVINGPLOTS) '_data.mat'];
        iter = importdata(fullfile(pname,fname));

        time = 1:iternum;
        juvenilephenos = iter(2:iternum+1,2);
        adultphenos = iter(2:iternum+1,4);
        juveniles = iter(2:iternum+1,1);
        adults = iter(2:iternum+1,3);

        for h = 1:iternum
            if any(juveniles{h}==0)
                ext  = find(juveniles{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end

            if any(adults{h}==0)
                ext  = find(adults{h}==0);
                juveniles{h}(ext)  = [];
                juvenilephenos{h}(ext)  = [];
                adults{h}(ext)  = [];
                adultphenos{h}(ext)  = [];
            end
        end

        % Time series of lineages scaled by population size
        figure(2)
        clf()
        set(gcf,'Position',[100 100 700 700])
        iter = iternum;
        for i=1:iternum
            x_axis = NaN(iternum,1);
            y_axis = NaN(iternum,1);
            x_axis(1:numel(juvenilephenos{i})) = juvenilephenos{i};
            y_axis(1:numel(adultphenos{i})) = adultphenos{i};
            z_axis = zeros(iternum,1);
            z_axis(1:numel(z_axis)) = i;
            sizes  = nan(iternum,1);
            sizes(1:numel(adults{i})) = sqrt(adults{i}/sum(adults{i}))*2e2;
            sizes = sqrt(sizes)*3e1;
            scatter3(x_axis,y_axis,z_axis,sizes,'MarkerFaceColor',col(ceil(i/iter*colres),:),'MarkerEdgeColor','k')
            hold on
        end
        %     ax=gca;
        %     ax.FontSize = 16;
        axis([-0.01 1 -0.01 1 0 iternum])
        xticks([0:0.5:1])
        yticks([0:0.5:1])
        zticks([0:100:iternum])
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',24)
        xlabel('\Theta_J')
        ylabel('\Theta_A')
        zlabel('Time')
        %legend('Lineages over time')
        sgtitle([' '])%Juvenile dominated. Time series of lineages scaled by population size. zj = ' num2str(zj) ', za = ' num2str(za)])
        view(0,90)
        pbaspect([1 1 1])
        figname = ['iter_' num2str(SAVINGPLOTS) '_lineages.fig'];
        saveas(gcf,fullfile(pname,figname))
        figname2 = ['iter_' num2str(SAVINGPLOTS) '_lineages.pdf'];
        print(fullfile(pname,figname2),'-bestfit','-dpdf')
    catch
    end
end

%%
ITERATIONNUMBER = 3; % CHANGE ACCORDING TO NUMBER OF ITERATIONS
pname = input('enter directory path name: ','s')
zj = -1.5;
za = -1.5;

for SAVINGPLOTS=1:ITERATIONNUMBER
    fname = ['iter_' num2str(SAVINGPLOTS) '_data.mat'];
    iter = importdata(fullfile(pname,fname));

    time = 1:iternum;
    juvenilephenos = iter(2:iternum+1,2);
    adultphenos = iter(2:iternum+1,4);
    juveniles = iter(2:iternum+1,1);
    adults = iter(2:iternum+1,3);

    for h = 1:iternum
        if any(juveniles{h}==0)
            ext  = find(juveniles{h}==0);
            juveniles{h}(ext)  = [];
            juvenilephenos{h}(ext)  = [];
            adults{h}(ext)  = [];
            adultphenos{h}(ext)  = [];
        end

        if any(adults{h}==0)
            ext  = find(adults{h}==0);
            juveniles{h}(ext)  = [];
            juvenilephenos{h}(ext)  = [];
            adults{h}(ext)  = [];
            adultphenos{h}(ext)  = [];
        end
    end

    % Time series of lineages scaled by population size
    figure(2)
    clf()
    set(gcf,'Position',[100 100 700 700])
    iter = iternum;
    for i=1:iternum
        x_axis = NaN(iternum,1);
        y_axis = NaN(iternum,1);
        x_axis(1:numel(juvenilephenos{i})) = juvenilephenos{i};
        y_axis(1:numel(adultphenos{i})) = adultphenos{i};
        z_axis = zeros(iternum,1);
        z_axis(1:numel(z_axis)) = i;
        sizes  = nan(iternum,1);
        sizes(1:numel(adults{i})) = sqrt(adults{i}/sum(adults{i}))*2e2;
        sizes = sqrt(sizes)*1e1;
        scatter3(x_axis,y_axis,z_axis,sizes,'MarkerFaceColor',col(ceil(i/iter*colres),:),'MarkerEdgeColor','k')
        hold on
    end
    axis([0 1 0 1 0 iternum])
    xticks([0:0.1:1])
    yticks([0:0.1:1])
    zticks([0:100:iternum])
    xlabel('\Theta_J')
    ylabel('\Theta_A')
    zlabel('Time')
    %legend('Lineages over time')
    sgtitle(['Juvenile dominated. Time series of lineages scaled by population size. zj = ' num2str(zj) ', za = ' num2str(za)])
    view(0,90)
    pbaspect([1 1 1])
    figname = ['iter_' num2str(SAVINGPLOTS) '_lineages.fig'];
    saveas(gcf,fullfile(pname,figname))
    figname2 = ['iter_' num2str(SAVINGPLOTS) '_lineages.pdf'];
    saveas(gcf,fullfile(pname,figname2),'pdf')
end

%% JUVENILES II
ITERATIONNUMBER = 1; % CHANGE ACCORDING TO NUMBER OF ITERATIONS
zj=-1; % CHANGE ACCORDING TO STRENGTH OF TRADE-OFF
za=-0.3; % CHANGE ACCORDING TO STRENGTH OF TRADE-OFF
inputdate = input('enter date in the folder name, format dd-Mon-yyyy: ','s')
for SAVINGPLOTS=1:ITERATIONNUMBER
    pname = '/Users/prm/Documents/MATLAB_scripts/Project2/Deterministic/isolated_mutations/Data';
    dname = ['JD_zj' num2str(zj) '_za' num2str(za) '_' inputdate '_sweep'];
    fname = ['iter_' num2str(SAVINGPLOTS) '_data.mat'];
    iter = importdata(fullfile(pname,dname,fname));

    xjparent = iter(2:iternum+1,1);
    xjparent = cell2mat(xjparent);
    xaparent = iter(2:iternum+1,2);
    xaparent = cell2mat(xaparent);
    invfitnessj = iter(2:iternum+1,3);
    invfitnessa = iter(2:iternum+1,4);
    for i = 2:iternum
        if invfitnessj{i} == []
            invfitnessj{i} = nan;
        else
            invfitnessj{i}=invfitnessj{i};%(1);
        end
        if invfitnessa{i} == []
            invfitnessa{i} = nan;
        else
            invfitnessa{i}=invfitnessa{i};%(1);
        end

        %invfitnessj{i}=invfitnessj{i}(1);
        %invfitnessa{i}=invfitnessa{i}(1);
    end
    invfitnessj = cell2mat(invfitnessj);
    invfitnessa = cell2mat(invfitnessa);
    xjm1 = iter(2:iternum+1,5);
    xjm1 = cell2mat(xjm1);
    xam1 = iter(2:iternum+1,6);
    xam1 = cell2mat(xam1);
    time = 1:iternum;
    juvenilephenos = iter(2:iternum+1,8);
    adultphenos = iter(2:iternum+1,10);
    juveniles = iter(2:iternum+1,7);
    adults = iter(2:iternum+1,9);

    % PLOT PARENT POP PHENOTYPE, MUTANT PHENOTYPE AND GRADIENT
    %     figure(1)
    %     clf()
    %     size = 60;
    %     subplot(2,1,1)
    %     scatter(time,xjparent,size,'MarkerEdgeColor','k','MarkerFaceColor',orange)
    %     hold on
    %     scatter(time,invfitnessj,size,'d','MarkerEdgeColor','k')
    %     hold on
    %     scatter(time,xjm1,size,'filled','s','MarkerEdgeColor','k','MarkerFaceColor',purple)
    %     xlabel('Time')
    %     hold on
    %     xticks([0:10:iternum])
    %     yticks([0:0.1:1])
    %     title(['Juvenile dominated. Time series of mutations and fitness gradient. zj = ' num2str(zj) ', za = ' num2str(za)])
    %     legend('\Theta_J parent pop','gradient J','\Theta_J mutant','Location','Best')
    %     subplot(2,1,2)
    %     scatter(time,xaparent,size,'MarkerEdgeColor','k','MarkerFaceColor',green)
    %     hold on
    %     scatter(time,invfitnessa,size,'d','MarkerEdgeColor','k')
    %     hold on
    %     scatter(time,xam1,size,'filled','s','MarkerEdgeColor','k','MarkerFaceColor',cyan)
    %     xlabel('Time')
    %     hold on
    %     xticks([0:10:iternum])
    %     yticks([0:0.1:1])
    %     legend('\Theta_A parent pop','gradient A','\Theta_A mutant','Location','Best')
    %     figname = ['iter_' num2str(SAVINGPLOTS) '_parent_mut_grad.fig'];
    %     saveas(gcf,fullfile(pname,dname,figname))
    %
    figure(2)
    clf()
    iter = iternum;
    for i=1:iternum
        x_axis = zeros(iternum,1);
        x_axis(1:numel(x_axis)) = i;
        y_axis = NaN(iternum,1);
        z_axis = NaN(iternum,1);
        y_axis(1:numel(juvenilephenos{i})) = juvenilephenos{i};
        z_axis(1:numel(adultphenos{i})) = adultphenos{i};
        sizes  = nan(iternum,1);
        sizes(1:numel(juveniles{i})) = juveniles{i};
        sizes = sqrt(sizes)*1e1;
        scatter3(x_axis,y_axis,z_axis,sizes,'MarkerFaceColor',col(ceil(i/iter*colres),:),'MarkerEdgeColor','k')
        hold on
    end
    axis([0 iternum 0 1 0 1])
    xticks([0:10:iternum])
    yticks([0:0.1:1])
    zticks([0:0.1:1])
    xlabel('Time')
    ylabel('\Theta_J')
    zlabel('\Theta_A')
    legend('Lineages over time')
    title(['Juvenile dominated. Time series of lineages scaled by population size. zj = ' num2str(zj) ', za = ' num2str(za)])
    figname = ['iter_' num2str(SAVINGPLOTS) '_lineages.fig'];
    saveas(gcf,fullfile(pname,dname,figname))

    %    Excel file
    iternum = length(invfitnessj);
    juvcol = nan(15,iternum);
    adcol  = nan(15,iternum);
    for i = 1:iternum
        jj = [xjparent(i); xjm1(i); invfitnessj(i); nan; juvenilephenos{i}];
        aa = [xaparent(i); xam1(i); invfitnessa(i); nan; adultphenos{i}];
        juvcol(1:numel(jj),i) = jj;
        adcol(1:numel(aa),i) = aa;
    end

    both = cell(1,iternum);
    for i = 1:iternum
        both{i} = [juvcol(:,i) adcol(:,i)];
    end

    % Order of columns: T_1 is juveniles and T_2 is adults
    % Order of rows: parent, mutant, R0, empty space, population vector
    T = table(both{1},'VariableNames',{'T1'});
    for i = 2:iternum
        Tnew = table(both{i},'VariableNames',{['T',num2str(i)]});
        T = [T, Tnew];
    end
    tabname = ['iter_' num2str(SAVINGPLOTS) '_mutation_sequence.xlsx'];
    writetable(T,fullfile(pname,dname,tabname))

end



%% EQUAL POPS
ITERATIONNUMBER = 10; % CHANGE ACCORDING TO NUMBER OF ITERATIONS
pname = input('enter directory path name: ','s')
zj = -0.5;
za = 0.5;

for SAVINGPLOTS=1:ITERATIONNUMBER
    fname = ['iter_' num2str(SAVINGPLOTS) '_data.mat'];
    iter = importdata(fullfile(pname,fname));

    time = 1:iternum;
    juvenilephenos = iter(2:iternum+1,2);
    adultphenos = iter(2:iternum+1,4);
    juveniles = iter(2:iternum+1,1);
    adults = iter(2:iternum+1,3);

    for h = 1:iternum
        if any(juveniles{h}==0)
            ext  = find(juveniles{h}==0);
            juveniles{h}(ext)  = [];
            juvenilephenos{h}(ext)  = [];
            adults{h}(ext)  = [];
            adultphenos{h}(ext)  = [];
        end

        if any(adults{h}==0)
            ext  = find(adults{h}==0);
            juveniles{h}(ext)  = [];
            juvenilephenos{h}(ext)  = [];
            adults{h}(ext)  = [];
            adultphenos{h}(ext)  = [];
        end
    end

    % Time series of lineages scaled by population size
    figure(2)
    clf()
    set(gcf,'Position',[100 100 700 700])
    iter = iternum;
    for i=1:iternum
        x_axis = NaN(iternum,1);
        y_axis = NaN(iternum,1);
        x_axis(1:numel(juvenilephenos{i})) = juvenilephenos{i};
        y_axis(1:numel(adultphenos{i})) = adultphenos{i};
        z_axis = zeros(iternum,1);
        z_axis(1:numel(z_axis)) = i;
        sizes  = nan(iternum,1);
        sizes(1:numel(adults{i})) = sqrt(adults{i}/sum(adults{i}))*2e2;
        sizes = sqrt(sizes)*1e1;
        scatter3(x_axis,y_axis,z_axis,sizes,'MarkerFaceColor',col(ceil(i/iter*colres),:),'MarkerEdgeColor','k')
        hold on
    end
    axis([-0.01 1 -0.01 1 0 iternum])
    xticks([0:0.5:1])
    yticks([0:0.5:1])
    zticks([0:100:iternum])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',24)
    xlabel('\Theta_J')
    ylabel('\Theta_A')
    zlabel('Time')
    %legend('Lineages over time')
    sgtitle([' '])% Time series of lineages scaled by population size. zj = ' num2str(zj) ', za = ' num2str(za)])
    view(0,90)
    pbaspect([1 1 1])
    figname = ['iter_' num2str(SAVINGPLOTS) '_lineages.fig'];
    saveas(gcf,fullfile(pname,figname))
    figname2 = ['iter_' num2str(SAVINGPLOTS) '_lineages.pdf'];
    print(fullfile(pname,figname2),'-bestfit','-dpdf')
end