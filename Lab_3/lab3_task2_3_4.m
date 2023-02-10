%Lab 3 Task 2
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%% x1-x2 grid to plot 2D filled contour plot 
x1_limits = [-10,5];
x2_limits = [-8,12];
n = 60;
[x1_plot,x2_plot] = meshgrid(linspace(x1_limits(1),x1_limits(2),n),linspace(x2_limits(1),x2_limits(2),n));

%%  comment
for i = 1:3 %comment
    %initial point, epsilon, gk, sk
    switch i 
        case 1
            x11 = 5;
            x21 = -5;
            epsilon = 0.01;
            gk = 0.5;
            sk = 5;
        case 2
            x11 = -5;
            x21 = 10;
            epsilon = 0.01;
            gk = 0.1;
            sk = 15;
        case 3
            x11 = 8;
            x21 = -10;
            epsilon = 0.01;
            gk = 0.2;
            sk = 0.1;
    end
    %run the algorithm
    [x1,x2] = StDes_proj_gk_const(x11,x21,epsilon,gk,sk,x1_limits,x2_limits);
    %plot the results
    figure();
    contourf(x1_plot,x2_plot,fx1x2(x1_plot,x2_plot),'--');
    title(['2D filled contour plot of f(x1,x2), Steepest Descent, gk = ',num2str(gk),', sk = ',num2str(sk),', epsilon = ',num2str(epsilon),', x1 = ',num2str(x11),' x2 = ',num2str(x21)]);
    hold on;
    plot(x1,x2,'red');
    scatter(x1,x2,'red');
    xlabel('x1-axis');
    ylabel('x2-axis');
    legend('f(x1,x2)','f(x1k,x2k)');
    hold off;
    figure()
    k = 1:length(x1);
    plot(k,fx1x2(x1,x2),'k');
    title(['f(x1,x2) values for different k, Steepest Descent, gk = ',num2str(gk),', sk = ',num2str(sk),', epsilon = ',num2str(epsilon),', x1 = ',num2str(x11),' x2 = ',num2str(x21)]);
    hold on;
    scatter(k,fx1x2(x1,x2),'red');
    xlabel('k');
    ylabel('f(x1k,x2k)');
    legend('f(x1,x2)','f(x1k,x2k)');
    hold off;
    if(i == 1) %extra plots for case 1
        figure()
        k = 1:length(x1);
        f_values = fx1x2(x1,x2);
        plot(k(1:length(k)/2),f_values(1:length(k)/2),'k');
        title(['f(x1,x2) values for different k, Steepest Descent, gk = ',num2str(gk),', sk = ',num2str(sk),', epsilon = ',num2str(epsilon),', x1 = ',num2str(x11),' x2 = ',num2str(x21)]);
        hold on;
        scatter(k(1:length(k)/2),f_values(1:length(k)/2),'red');
        xlabel('k');
        ylabel('f(x1k,x2k)');
        legend('f(x1,x2)','f(x1k,x2k)');
        hold off;
        figure()
        plot(k(length(k)/2+1:length(k)),f_values(length(k)/2+1:length(k)),'k');
        title(['f(x1,x2) values for different k, Steepest Descent, gk = ',num2str(gk),', sk = ',num2str(sk),', epsilon = ',num2str(epsilon),', x1 = ',num2str(x11),' x2 = ',num2str(x21)]);
        hold on;
        scatter(k(length(k)/2+1:length(k)),f_values(length(k)/2+1:length(k)),'red');
        xlabel('k');
        ylabel('f(x1k,x2k)');
        legend('f(x1,x2)','f(x1k,x2k)');
        hold off;
    end
end

%% function to execute Steepest Descent method with projection and gk constant
function [x1,x2] = StDes_proj_gk_const(x11,x21,epsilon,gk,sk,x1_limits,x2_limits)
    if(x11<x1_limits(1))
       fprintf('\nInitial point is outside of the limits - x1 = %i\n',x11);
       x11 = x1_limits(1);
    end
    if(x11>x1_limits(2))
        fprintf('\nInitial point is outside of the limits - x1 = %i\n',x11);
        x11 = x1_limits(2);
    end
    if(x21<x2_limits(1))
        fprintf('\nInitial point is outside of the limits - x2 = %i\n',x21);
        x21 = x2_limits(1);
    end
    if(x21>x2_limits(2))
        fprintf('\nInitial point is outside of the limits - x2 = %i\n',x21);
        x21 = x2_limits(2);
    end
    x1(1) = x11;
    x2(1) = x21;
    k = 1;
    while norm(grad_fx1x2(x1(k),x2(k))) > epsilon 
        d(:,k) = -grad_fx1x2(x1(k),x2(k));
        x1k_bar = x1(k)+sk*d(1,k);
        x2k_bar = x2(k)+sk*d(2,k);
        if(x1k_bar<x1_limits(1))
            x1k_bar = x1_limits(1);
        end
        if(x1k_bar>x1_limits(2))
            x1k_bar = x1_limits(2);
        end
        if(x2k_bar<x2_limits(1))
            x2k_bar = x2_limits(1);
        end
        if(x2k_bar>x2_limits(2))
            x2k_bar = x2_limits(2);
        end
        x1(k+1) = x1(k)+gk*(x1k_bar - x1(k));
        x2(k+1) = x2(k)+gk*(x2k_bar - x2(k));
        k = k + 1;
        if(k==500)
            break;
        end  
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