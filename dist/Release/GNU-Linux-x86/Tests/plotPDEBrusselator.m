clear;
run(['PDEBrusselator/solution.m']);

n = size(y,2);
n = n/2;

u = y(:,1:n);
v = y(:,(n+1):end);
x = linspace(0,1,n+2);
x = x(2:(end-1));

surf(x,t,u);
shading interp;
% figure;
% surf(x,t,v);

