clear;

% testname = 'NeuronCable';
% testname = 'OneDimensionalTest1';
% testname = 'TwoDimensionalBrusselator';
% testname = 'PDEBrusselator';
testname = 'Krogh10';

solvername = 'DROCK2';

tol_pow = [2, 3, 4, 5, 6, 7, 8];
dt_pow = [1, 2, 4, 6, 8, 10, 12];
dt = 1./2.^dt_pow;

refsol = csvread([testname '/reference.txt']);
neqn = numel(refsol);
err_adap = zeros(numel(tol_pow),1);
cost_adap = zeros(numel(tol_pow),1);
for i=1:numel(tol_pow)
    err_adap(i) = max(abs(csvread([testname '/EfficiencyTest_' solvername '_tol_10e_' int2str(tol_pow(i)) '.txt'])-refsol));
    tmp = csvread([testname '/EfficiencyTest_' solvername '_tol_10e_' int2str(tol_pow(i)) '_statistics.txt']);
    cost_adap(i) = tmp(1); % 1 for time, 2 for function evaluations
end

err_fix = zeros(numel(dt_pow),1);
cost_fix = zeros(numel(dt_pow),1);
for i=1:numel(tol_pow)
    err_fix(i) = max(abs(csvread([testname '/EfficiencyTest_' solvername '_dt_2e_' int2str(dt_pow(i)) '.txt'])-refsol));
    tmp = csvread([testname '/EfficiencyTest_' solvername '_dt_2e_' int2str(dt_pow(i)) '_statistics.txt']);
    cost_fix(i) = tmp(1); % 1 for time, 2 for function evaluations
end

figure;
loglog(err_adap, cost_adap,'b');
hold on;
loglog(err_fix,cost_fix,'r');