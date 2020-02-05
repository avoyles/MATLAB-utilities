function [X_norm, Y_norm] = real2norm(x_plot, y_plot)

axPos = get(gca,'Position'); %# gca gets the handle to the current axes
xMinMax = xlim;
yMinMax = ylim;
X_norm = axPos(1) + ((x_plot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
Y_norm = axPos(2) + ((y_plot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);

end