function nummarkers(h,num)

% NUMMARKERS takes a vector of line handles in h
% and reduces the number of plot markers on the lines
% to num. This is useful for closely sampled data.
%
% example:
% t = 0:0.01:pi;
% p = plot(t,sin(t),'-*',t,cos(t),'r-o');
% nummarkers(p,10);
% legend('sin(t)','cos(t)')
%

% Magnus Sundberg Feb 08, 2001

for n = 1:length(h)
    if strcmp(get(h(n),'type'),'line')
        axes(get(h(n),'parent'));
        x = get(h(n),'xdata');
        y = get(h(n),'ydata');
        t = 1:length(x);
        s = [0 cumsum(sqrt(diff(x).^2+diff(y).^2))];
        si = (0:num-1)*s(end)/(num-1);
        ti = round(interp1(s,t,si));
        xi = x(ti);
        yi = y(ti);
        marker = get(h(n),'marker');
        color = get(h(n),'color');
        style = get(h(n),'linestyle');
        % make a line with just the markers
        set(line(xi,yi),'marker',marker,'linestyle','none','color',color);
        % make a copy of the old line with no markers
        set(line(x,y),'marker','none','linestyle',style,'color',color);
        % set the x- and ydata of the old line to [], this tricks legend to keep on working
        set(h(n),'xdata',[],'ydata',[]);
    end
end