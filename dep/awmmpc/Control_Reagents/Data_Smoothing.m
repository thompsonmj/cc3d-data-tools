%% Generate new data structure
%{
clear;clc;
load('data.mat');
data(7).y_obs = data(1).y_obs;
data(7).t_obs = data(1).t_obs;
data(13).y_obs = data(1).y_obs;
data(13).t_obs = data(1).t_obs;
data_temp = data;
Nexpt = numel(data);
ti = 0:1:30;
p = 0.6;

fig = figure; set(fig,'color','w');
for i = 1:Nexpt
    subplot(round(sqrt(Nexpt)),ceil(sqrt(Nexpt)),i);
    N = size(data(i).y_obs{1},1);
    y = data(i).y_obs{1};
    t = data(i).t_obs;
%     ti = t;
    yi = csaps(t,mean(y,1),p,ti);
    plot(t,y,'ro',ti,yi,'b-');
    data_temp(i).y_obs{1} = yi;
    data_temp(i).t_obs = ti;
end
data = data_temp;
% save('data_smooth.mat','data');
%}

%% Plot 2D projections
% %{
clear;clc;close all;

fontname = 'helvetica'; fontsize = 8;
linewidth = 1; markersize = 7;
load('data_smooth.mat');
data_smooth = data;
load('data.mat');
data(7) = data(1); data(13) = data(1);
D{1} = data(1:6); D{2} = data(7:12); D{3} = data(13:16);
DS{1} = data_smooth(1:6); DS{2} = data_smooth(7:12); DS{3} = data_smooth(13:16);
co = {'y','b','g','r','c','m','y'};
% co = {'b','b','b','b','b','b'};
sh = {'v','s','x','o','+','^'};
c = [50,10];
figdim = [3,2.25];

figdim(1) = max(123,8+96*figdim(1));
figdim(2) = 82+96*figdim(2);
for j = 1:3
    fig = figure(j); clf; set(fig,'color','w'); hold on;
%     fig = figure(1); set(fig,'color','w'); hold on;
    rect = get(fig,'outerposition'); set(fig,'outerposition',[rect(1),rect(2)/1.25,figdim]);
    set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel('Time (minutes)'); ylabel('Normalized [pErk]');
    data = D{j}; data_smooth = DS{j}; Nexpt = numel(data);
    h = []; txt = []; set(gca,'fontname',fontname,'fontsize',fontsize);%16);
    for i = 1:Nexpt
        u = data(i).u_admin;
        t = data(i).t_obs;
        y = data(i).y_obs{1};
        ymean = mean(y,1);
        E = std(y,[],1)/sqrt(size(y,1));
        yi = data_smooth(i).y_obs{1};% yi(4:6) = [0.8,0.92,0.86];
        ti = data_smooth(i).t_obs;
        plot(ti,yi,[co{i},'-'],'linewidth',linewidth);
        h1 = errorbar(t,ymean,E,[co{i},sh{i}],'linewidth',linewidth,'markersize',markersize);
        if ismember(j,[1,2]), u_lgnd = num2str(c(j)*u(j)); ...
        else u_lgnd = ['[',num2str(c(1)*u(1)),', ',num2str(c(2)*u(2)),']']; end
        [~,~,h,txt] = legend([h;h1],[txt,[u_lgnd,' \muM']],'fontsize',fontsize);%12);
    end
    axis([0,30,0,1.400000000001]);% box on; grid on;
end
%}

%% Plot 3D
% %{
clear;clc;%close all;

fontname = 'times new roman'; fontsize = 12;
linewidth = 3; markersize = 10;
load('data_smooth.mat');
data_smooth = data;
load('data.mat');
data(7) = data(1);
D{1} = data(1:6); D{2} = data(7:12);
DS{1} = data_smooth(1:6); DS{2} = data_smooth(7:12);
co = {'b','g','r','c','m','y'}; c = [50,10];
for j = 1:2
    fig = figure(j); clf; set(fig,'color','w');
    set(gca,'fontname',fontname,'fontsize',fontsize);
    data = D{j}; data_smooth = DS{j}; Nexpt = numel(data);
    h = []; txt = [];
    for i = 1:Nexpt
        u = data(i).u_admin;
        t = data(i).t_obs;
        y = data(i).y_obs{1};
        ymean = mean(y,1);
        E = std(y,[],1)/sqrt(size(y,1));
        yi = data_smooth(i).y_obs{1};
        ti = data_smooth(i).t_obs;
        h1 = plot3(ti,c(j)*u(j)*ones(size(ti)),yi,[co{i},'-'],'linewidth',linewidth); hold on;
        plot3(t,c(j)*u(j)*ones(size(t)),ymean,[co{i},'o'],'linewidth',linewidth,'markersize',markersize);
        z = [ymean+E;ymean-E];
        for k = 1:numel(t)
            plot3([t(k);t(k)],c(j)*u(j)*ones(2,1),z(:,k),[co{i},'-'],'linewidth',linewidth);
        end
        u_lgnd = num2str(c(j)*u(j));
        [~,~,h,txt] = legend([h;h1],[txt,[u_lgnd,' \muM']]);
    end
    xlabel('Time (minutes)'); ylabel('Input Concentration (\muM)');
    zlabel('Normalized [pErk]'); axis([0,30,0,c(j),0,1.4]); grid on;
end
%}