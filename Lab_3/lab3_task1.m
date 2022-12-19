%Lab 3 Task 1
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%% x1-x2 grid to plot 2D filled contour plot 
start = -12;
stop = 12;
n = 60;
[x1_plot,x2_plot] = meshgrid(linspace(start,stop,n),linspace(start,stop,n));

%% initial point
x11 = -10;
x21 = -8;

%% Running Steepest Descent method for initial point (x11,x21) and for every gk case
for i = 1:4
    %testing different gk values / 4 cases : 0.1, 0.3, 3, 5
    switch i
        case 1
            gk = 0.1; 
        case 2
            gk = 0.3; 
        case 3
            gk = 3; 
        case 4
            gk = 5; 
    end
    epsilon = 0.001;
    [x1,x2] = StDes_gk_const(x11,x21,gk,epsilon);
    %plot the results
    figure();
    contourf(x1_plot,x2_plot,fx1x2(x1_plot,x2_plot),'--');
    title(['2D filled contour plot of f(x1,x2), Steepest Descent, gk = ',num2str(gk),', x1 = ',num2str(x11),' x2 = ',num2str(x21)]);
    hold on;
    plot(x1,x2,'red');
    scatter(x1,x2,'red');
    xlabel('x1 axis');
    ylabel('x2 axis');
    legend('f(x1,x2)','f(x1k,x2k)');
    hold off;
    figure()
    k = 1:length(x1);
    plot(k,fx1x2(x1,x2),'k');
    title(['f(x1,x2) values for different k, Steepest Descent, gk = ',num2str(gk),', x1 = ',num2str(x11),' x2 = ',num2str(x21)]);
    hold on;
    scatter(k,fx1x2(x1,x2),'red');
    xlabel('k');
    ylabel('f(x1k,x2k)');
    legend('f(x1,x2)','f(x1k,x2k)');
    hold off;
end

%% function to execute Steepest Descent method with gk constant
function [x1,x2] = StDes_gk_const(x11,x21,gk,epsilon)
    x1(1) = x11;
    x2(1) = x21;
    k = 1;
    while norm(grad_fx1x2(x1(k),x2(k))) > epsilon
        d(:,k) = -grad_fx1x2(x1(k),x2(k));
        x1(k+1) = x1(k)+gk*d(1,k);
        x2(k+1) = x2(k)+gk*d(2,k);
        k = k + 1;       
    end
end

%% function to calculate the f(x1,x2) at a point (fx1k,x2k)
function fx1kx2k = fx1x2(x1k,x2k)
    syms x1 x2;
    fx1x2 = 1/3*x1.^2 + 3*x2.^2;
    fx1kx2k = subs(fx1x2,{x1,x2},{x1k,x2k});
end
%% function to calculate the gradient of f(x1,x2) at X1,X2
function [grad_fX1X2] = grad_fx1x2(X1,X2)
    syms x1 x2;
    fx1x2 = 1/3*x1.^2 + 3*x2.^2;
    grad_fx1x2 = gradient(fx1x2,[x1,x2]);
    grad_fX1X2 = subs(grad_fx1x2,[x1 x2],{X1,X2});
end



