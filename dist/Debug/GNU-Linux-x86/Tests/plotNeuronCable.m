rehash ;
clear;
pack;
clc;
% close all;


run(['NeuronCable/solution.m']);

nsol=numel(t);

figure;
hold on;
for i=1:nsol
    plot(y(i,:));
end


% y = csvread(['NeuronCable/solution.txt']);
% plot(y);


