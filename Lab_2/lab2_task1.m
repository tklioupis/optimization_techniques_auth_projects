%Lab 2 Task 1
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%% plotting f(x,y) = x^5*exp(-x^2-y^2)
start = -5;
stop = 5;
n = 60;
[x,y] = meshgrid(linspace(start,stop,n),linspace(start,stop,n));
fxy = x.^5.*exp(-x.^2-y.^2);
surf(x,y,fxy);
title('f(x,y) = x^5 * exp(-x^2-y^2)');
view(-15,15);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('f(x,y)');
