clear;
testname = 'SdeOneDimensionalTest2';
% testname = 'SdeLinear';

solvername = 'SROCK2';
% solvername = 'MT';

p = [0,1,2,3,4];
dt = 1./2.^p;

err = zeros(numel(p),8);
for i=1:numel(p)
    err(i,:) = csvread([testname '/TimeConvTest_' solvername '_dt_' int2str(p(i)) '_mt_error.txt']);
end

loglog(dt,err(:,1));
hold on;
loglog(dt,err(:,1),'.');
loglog(dt,1e-1*dt.^2,'--k');
loglog(dt,dt,'k');
