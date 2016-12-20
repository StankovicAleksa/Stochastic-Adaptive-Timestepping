sol_1e_6;
y6=y(:,2);
sol_1e_5;
y5=y(:,2);
sol_1e_4;
y4=y(:,2);
sol_1e_3;
y3=y(:,2);
sol_1e_2;
y2=y(:,2);
sol_1e_1;
y1=y(:,2);
sol_1e_0;
y0=y(:,2);

err(1)=max(abs(y6-y0));
err(2)=max(abs(y6-y1));
err(3)=max(abs(y6-y2));
err(4)=max(abs(y6-y3));
err(5)=max(abs(y6-y4));
err(6)=max(abs(y6-y5));
dt = [1e-0 1e-1 1e-2 1e-3 1e-4 1e-5];

loglog(dt,err);
hold on;
loglog(dt,dt.^2,'--k');