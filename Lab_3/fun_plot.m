%script to plot the function of Lab3
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%% plotting f(x1,x2) = 1/3*x1^2 + 3*x2^2
start = -12;
stop = 12;
n = 60;
[x1,x2] = meshgrid(linspace(start,stop,n),linspace(start,stop,n));
f = 1/3*x1.^2 + 3*x2.^2;
surf(x1,x2,f);
title('f(x1,x2) = 1/3*x1^2 + 3*x2^2');
view(-30,15);
xlabel('x1-axis');
ylabel('x2-axis');
zlabel('f(x1,x2)');
