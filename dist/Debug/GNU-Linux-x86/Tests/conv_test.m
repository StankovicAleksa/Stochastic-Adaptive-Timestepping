
testname = 'NeuronCable';
testname = 'OdeOneDimensionalTest1';
testname = 'Krogh10';
testname = 'PDEBrusselator';
testname = 'TwoDimensionalBrusselator';

solvername = 'ROCK2';

p = [2,4,6,8,10,12,14,16];
p = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
dt = 1./2.^p;

refsol = csvread([testname '/reference.txt']);
neqn = numel(refsol);
err = zeros(numel(p),1);
for i=1:numel(p)
    err(i) = max(abs(csvread([testname '/TimeConvTest_' solvername '_dt_' int2str(p(i)) '.txt'])-refsol));
end

loglog(dt,err);
hold on;
loglog(dt,dt.^2,'--k');