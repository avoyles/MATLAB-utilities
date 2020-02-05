function x = horz(y,xf,l,x0);
switch nargin
    case 2
        x = linspace(0,xf,1000);
    case 3
        x = linspace(0,xf,l);
    case 4
        x = linspace(x0,xf,l);
end