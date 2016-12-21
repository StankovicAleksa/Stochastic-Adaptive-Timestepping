clear;
run(['TwoDimensionalBrusselator/solution.m']);

plot(y(:,1),y(:,2));

figure;
plot(y(:,1));
hold on;
plot(y(:,2));

