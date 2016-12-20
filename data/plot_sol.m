clear;
mesh;
sol;

nsol=numel(t);

figure;
hold on;
for i=1:nsol
    plot(x,y(:,i));
end