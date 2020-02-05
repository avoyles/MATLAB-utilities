function [X Y] = ellipse(Xc,Yc,dX,dY,N_p)
%calculates an elipse

theta=linspace(0,2*pi,N_p);


X=Xc+dX*cos(theta)-dY*sin(theta);
Y=Yc+dX*cos(theta)+dY*sin(theta);

X=Xc+dX*cos(theta);
Y=Yc+dY*sin(theta);

% for i=1:length(Xc);
%     
%     X(i)=Xc(i)+dX(i)*cos(theta)-dY(i)*sin(theta);
%     Y(i)=Yc(i)+dX(i)*cos(theta)+dY(i)*sin(theta);
% 
%     X(i)=Xc(i)+dX(i)*cos(theta);
%     Y(i)=Yc(i)+dY(i)*sin(theta);
% end