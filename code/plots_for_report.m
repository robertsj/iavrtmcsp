% loads and plots some data
function plots_for_report
clc
xx = 0.25:0.5:10;


% printing for three group source detector
load slab_3g_analog.mat 
  disp('---------------------')
  out.phi2(20), out.re2(20), out.fom2(20), out.t,  2*snout.phiFavg(20)
load slab_3g_impcapt.mat
  disp('---------------------')
  out.phi2(20), out.re2(20), out.fom2(20), out.t,  2*snout.phiFavg(20)
load slab_3g_geomsplt.mat
  disp('---------------------')
  out.phi2(20), out.re2(20), out.fom2(20), out.t,  2*snout.phiFavg(20)
load slab_3g_stcadis_cell_20.mat
  disp('---------------------')
  out.phi2(20), out.re2(20), out.fom2(20), out.t,  2*snout.phiFavg(20)
load slab_3g_coopers.mat
  disp('---------------------')
  out.phi2(20), out.re2(20), out.fom2(20), out.t,  2*snout.phiFavg(20)
return




% one-group global statistics
figure(1)
title('Track Length Estimated Flux -- Relative Error')
xlabel('x [cm]'), ylabel('relative error')
hold on
for i = 1:3
    if i == 1
        load slab_1g_fwcadis_fullB.mat
    elseif i == 2
        load slab_1g_coopersB.mat
    else
        load slab_1g_stcadis_cell_20B.mat
    end
    if i==1
        plot(xx,out.re2,'-','Color',plotcolor(1),'LineWidth',2)
    elseif i == 2
        plot(xx,out.re2,'--','Color',plotcolor(2),'LineWidth',2)
    else
        plot(xx,out.re2,'-.','Color',plotcolor(3),'LineWidth',2)
    end
end
legend('FW-CADIS','Cooper','CADIS',0)
grid on

% one-group Monte Carlo density
figure(2)
title('Monte Carlo Particle Density')
xlabel('x [cm]'), ylabel('m(x,E_g) [mcp/cm^3-sp]')
hold on
for i = 1:3
    if i == 1
        load slab_1g_fwcadis_fullB.mat
    elseif i == 2
        load slab_1g_coopersB.mat
    else
        load slab_1g_stcadis_cell_20B.mat
    end
    for g = 1:input.numg
        fig = errorbar(xx,out.phiMC(:,g)/2,...
            out.reMC(:,g).*out.phiMC(:,g)/2,mark(i),'LineWidth',2);
        set(fig,'Color',plotcolor(i),'MarkerEdgeColor','k',...
            'MarkerFaceColor',plotcolor(i),...
            'MarkerSize',6);
    end
end
legend('FW-CADIS','Cooper','CADIS',0)
grid on

% three-group global statistics
figure(3)
title('Track Length Estimated Flux -- Relative Error')
xlabel('x [cm]'), ylabel('relative error')
hold on
for g = 1:1:3
    for i = 1:3
        if i == 1
            load slab_3g_fwcadis_full.mat
        elseif i == 2
            load slab_3g_coopers.mat
        else
            load slab_3g_stcadis_cell_20.mat
        end
        if i == 1
            plot(xx,out.re2(:,g),mark(g),'Color',plotcolor(i),'LineWidth',2)
            lab(1,:)=(['FW-CADIS']);
        elseif i == 2
            plot(xx,out.re2(:,g),mark(g),'Color',plotcolor(i),'LineWidth',2)
            lab(2,:)=(['Cooper  ']);
        else
            plot(xx,out.re2(:,g),mark(g),'Color',plotcolor(i),'LineWidth',2)
            lab(3,:)=(['CADIS   ']);
        end
    end
end
legend(lab,0)
grid on

% three-group Monte Carlo density
figure(4)
title('Monte Carlo Particle Density')
xlabel('x [cm]'), ylabel('m(x,E_g) [mcp/cm^3-sp]')
hold on
for g = 1:3
    for i = 1:3
        if i == 1
            load slab_3g_fwcadis_full.mat
        elseif i == 2
            load slab_3g_coopers.mat
        else
            load slab_3g_stcadis_cell_20.mat
        end
        fig = errorbar(xx,out.phiMC(:,g)/2,...
            out.reMC(:,g).*out.phiMC(:,g)/2,mark(g),'LineWidth',2);
        set(fig,'Color',plotcolor(i),'MarkerEdgeColor','k',...
            'MarkerFaceColor',plotcolor(i),...
            'MarkerSize',6);
    end
end
legend(lab,0)
grid on



end

function lala = mark(i)
    switch i
        case 1
            lala = '-'; % blue
        case 2
            lala = '--'; % nice green
        case 3
            lala = '-.'; % red
    end
end


function color = plotcolor(g)
    % this function is a hard-coded color map for the different
    % flux groups.  I've accounted for up to 8 groups.
    % set(ur,'Color',[1 0.7 0.2],'LineWidth',2);
    switch g
        case 1
            color = [0.0 0.0 1.0]; % blue
        case 2
            color = [0.0 0.8 0.2]; % nice green
        case 3
            color = [1.0 0.0 0.0]; % red
        case 4
            color = [0.4 0.0 0.6]; % purple
        case 5
            color = [0.9 0.4 0.0]; % orange
        case 6
            color = [0.5 0.2 0.0]; % brown
        case 7
            color = [0.0 0.8 0.6]; % turquoise
        case 8
            color = [0.7 0.6 0.0]; % gold
        otherwise
            
    end
end