function [x_best,y_best,cost] = getInputs(nncOut,w)

x = cell2mat(struct2cell([nncOut.xVals{:}]')');
y = [nncOut.yVals{:}]';
x = x(nncOut.gblPrto,:);
y = y(nncOut.gblPrto,:);
one = ones(size(y,1),1);
yUtopia = nncOut.yUtopia(:)';
yNadir = nncOut.yNadir(:)';
L = yNadir - yUtopia; L(~L) = 1;
yn = (y - one*yUtopia)./(one*L);

% Norm of vector weighted by Akaike weights
w = w/sum(w); w = one*w(:)';
ynL2norm = sum(w.*(yn.^2),2);
[~,minInd] = min(ynL2norm);

x_best = x(minInd,:);
y_best = y(minInd,:);
cost = ynL2norm;
%{
w_old = w(1,:)
x_best
y_best
yn_min = yn(minInd,:)
ynL2norm_min = min(ynL2norm)
figure(1); clf; axis([0,1,0,1]); hold on; caxis([0,1/epsilon]);
scatter(x(:,1),x(:,2),40,ynL2norm,'filled');
plot(x_best(1),x_best(2),'m+','linewidth',3,'markersize',15);
xlabel('u_1'); ylabel('u_2');
figure(2); clf; hold on;
scatter3(yn(:,1),yn(:,2),yn(:,3),40,ynL2norm,'filled');
plot3(yn(minInd,1),yn(minInd,2),yn(minInd,3),'m+','linewidth',3,'markersize',15);
xlabel('Objective 1'); ylabel('Objective 2'); zlabel('Objective 3');
colorbar
%}
end