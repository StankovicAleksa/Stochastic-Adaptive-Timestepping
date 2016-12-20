clear;

mtdata_1e_0;
means(1,:) = mean;
dts(1)= dt;
mtdata_1e_1;
means(2,:) = mean;
dts(2)= dt;
mtdata_1e_2;
means(3,:) = mean;
dts(3)= dt;
mtdata_1e_3;
means(4,:) = mean;
dts(4)= dt;
mtdata_1e_4;
means(5,:) = mean;
dts(5)= dt;
mtdata_1e_5;
means(6,:) = mean;

err = zeros(size(means,1)-1,size(means,2));
for i=1:5
    err(i,:) = abs(means(5,:)-means(i,:));
end

subplot(2,2,1);
loglog(dts,err(:,1));
hold on;
loglog(dts,dts.^2,'--k');
leg=legend('$\int y dx$','$\mathcal{O}(\Delta t^2)$');
set(leg,'interpreter','latex');

subplot(2,2,2);
loglog(dts,err(:,2));
hold on;
loglog(dts,dts.^2,'--k');
leg=legend('$\int \sin(y) dx$','$\mathcal{O}(\Delta t^2)$');
set(leg,'interpreter','latex');

subplot(2,2,3);
loglog(dts,err(:,3));
hold on;
loglog(dts,dts.^2,'--k');
leg=legend('$\int y^2 dx$','$\mathcal{O}(\Delta t^2)$');
set(leg,'interpreter','latex');

subplot(2,2,4);
loglog(dts,err(:,4));
hold on;
loglog(dts,dts.^2,'--k');
leg=legend('$\int e^y dx$','$\mathcal{O}(\Delta t^2)$');
set(leg,'interpreter','latex');
