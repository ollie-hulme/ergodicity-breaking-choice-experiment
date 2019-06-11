function h = raincloud_plot(X,cl)

[a,b] = ksdensity(X);

wdth = 0.8; % width of boxplot
% TODO, should probably be some percentage of max.height of kernel density plot

% density plot
h{1} = area(b,a); hold on
set(h{1}, 'FaceColor', cl);
set(h{1}, 'EdgeColor', [0.1 0.1 0.1]);
set(h{1}, 'LineWidth', 1);

% make some space under the density plot for the boxplot
yl = get(gca,'YLim');
set(gca,'YLim',[-2 yl(2)]);

% jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% info for making boxplot
Y = quantile(X,[0.25 0.75 0.5 0.02 0.98]);

% 'box' of 'boxplot'
h{2} = rectangle('Position',[Y(1) -1-(wdth*0.5) Y(2)-Y(1) wdth]);
set(h{2},'EdgeColor','k')
set(h{2},'LineWidth',1);
% could also set 'FaceColor' here as Micah does, but I prefer without

% mean line
h{3} = line([Y(3) Y(3)],[-1.2 -0.8],'col','k','LineWidth',1);

% whiskers
h{4} = line([Y(2) Y(5)], [-1 -1],'col','k','LineWidth',1);
h{5} = line([Y(1) Y(4)],[-1 -1],'col','k','LineWidth',1);

% raindrops
h{3} = scatter(X,jit - 1);
h{3}.SizeData = 5;
h{3}.MarkerFaceColor = cl;
h{3}.MarkerEdgeColor = 'none';



