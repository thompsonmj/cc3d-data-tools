%%
clear;clc;%close all;

t = 0:0.1:3000;
[x0,P] = TCell_Lipn;% x0 = P.x0;
x0(1) = 0;
% x0([31,33,34,36,37]) = [0,0,0,0,0];
% load('P_opt.mat','P_L'); P = P_L;
% out = [36,37]; C = [1,1];
out = 1:length(x0); C = eye(numel(out));
Ny = size(C,1); %P.m2 = 0.0469;
% outLabel = {'CD45','ZAP70','MEK','ERK'};
[T,X] = ode15s(@TCell_Lipn,t,x0,[],P,[],[]);
Y = C*X(:,out)';

fig = figure; set(fig,'color','w');
for j = 1:Ny
    subplot(round(sqrt(Ny)),ceil(sqrt(Ny)),j);
    hold on; plot(T,Y(j,:),'r','linewidth',2);% ylabel(outLabel{j});
    ylabel(num2str(j));
end
xeq = X(end,:)';
figure; hold on;
plot(log10(x0+1),'bo','linewidth',2)
plot(log10(xeq+1),'rx','linewidth',2)
dx = TCell_Lipn(0,x0,P,[],[])
norm(dx,2)
dx = TCell_Lipn(0,xeq,P,[],[])
norm(dx,2)


%%
xeq(1) = 1e4;
x0(1) = 1e4;
t = 0:0.1:30;
out = [36,37]; C = [1,1]; Ny = size(C,1);
[T1,X1] = ode15s(@TCell_Lipn,t,xeq,[],P,[],[]);
[T2,X2] = ode15s(@TCell_Lipn,t,x0,[],P,[],[]);
Y = C*X2(:,out)';
fig = figure; set(fig,'color','w');
for j = 1:Ny
    subplot(round(sqrt(Ny)),ceil(sqrt(Ny)),j);
    hold on;
%     plot(T2,Y(j,:),'r','linewidth',2);
    plot(T1,X1(:,out(j)),'r','linewidth',2);
    plot(T2,X2(:,out(j)),'k','linewidth',2);
%     plotyy(T1,X1(:,out(j)),T2,X2(:,out(j)));
%     ylabel(outLabel{j});
    ylabel(num2str(j));
end