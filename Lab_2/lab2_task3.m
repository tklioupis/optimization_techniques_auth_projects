%Lab 2 Task 3
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%% x-y grid to plot 2D filled contour plot 
start = -3;
stop = 3;
n = 60;
[x_plot,y_plot] = meshgrid(linspace(start,stop,n),linspace(start,stop,n));

%% Running Newton method for every initial point (x1,y1) and for every gk case
for i = 1:3
    %initialization of (x1,y1) / 3 cases : (0,0) (-1,1) (1,-1)
    switch i
        case 1
            x1 = 0;
            y1 = 0;
        case 2
            x1 = -1;
            y1 = 1;
        case 3
            x1 = 1;
            y1 = -1;
    end
    epsilon = 0.001;
    for j = 1:3 %3 cases for gk: const, minimizes f(xk+gk*dk), armijo
        switch j
            case 1 %gk const
                gk = 0.1;
                [x,y] = newton_gk_const(x1,y1,gk,epsilon);
                %plot the results
                figure();
                contourf(x_plot,y_plot,fxy(x_plot,y_plot),'--');
                title(['2D filled contour plot of f(x,y), Newton, gk const,',' x1 = ',num2str(x1),' y1 = ',num2str(y1)]);
                hold on;
                plot(x,y,'red');
                scatter(x,y,'red');
                xlabel('x axis');
                ylabel('y axis');
                legend('f(x,y)','f(xk,yk)');
                hold off;
                figure()
                k = 1:length(x);
                plot(k,fxy(x,y),'k');
                title(['f(x,y) values for different k, Newton, gk const,',' x1 = ',num2str(x1),' y1 = ',num2str(y1)]);
                hold on;
                scatter(k,fxy(x,y),'red');
                xlabel('k');
                ylabel('f(xk,yk)');
                legend('f(x,y)','f(xk,yk)');
                hold off;
            case 2 %gk that minimizes f(xk+gk*dk)
                [x,y] = newton_gk_minf(x1,y1,epsilon);
                %plot the results
                figure();
                contourf(x_plot,y_plot,fxy(x_plot,y_plot),'--');
                title(['2D filled contour plot of f(x,y), Newton, gk that minimizes f(xk+gk*dk),',' x1 = ',num2str(x1),' y1 = ',num2str(y1)]);
                hold on;
                plot(x,y,'red');
                scatter(x,y,'red');
                xlabel('x axis');
                ylabel('y axis');
                legend('f(x,y)','f(xk,yk)');
                hold off;
                figure()
                k = 1:length(x);
                plot(k,fxy(x,y),'k');
                title(['f(x,y) values for different k, Newton, gk that minimizes f(xk+gk*dk),',' x1 = ',num2str(x1),' y1 = ',num2str(y1)]);
                hold on;
                scatter(k,fxy(x,y),'red');
                xlabel('k');
                ylabel('f(xk,yk)');
                legend('f(x,y)','f(xk,yk)');
                hold off;
            case 3 %gk armijo
                [x,y] = newton_gk_armijo(x1,y1,epsilon);
                %plot the results
                figure();
                contourf(x_plot,y_plot,fxy(x_plot,y_plot),'--');
                title(['2D filled contour plot of f(x,y), Newton, gk armijo,',' x1 = ',num2str(x1),' y1 = ',num2str(y1)]);
                hold on;
                plot(x,y,'red');
                scatter(x,y,'red');
                xlabel('x axis');
                ylabel('y axis');
                legend('f(x,y)','f(xk,yk)');
                hold off;
                figure()
                k = 1:length(x);
                plot(k,fxy(x,y),'k');
                title(['f(x,y) values for different k, Newton, gk armijo,',' x1 = ',num2str(x1),' y1 = ',num2str(y1)]);
                hold on;
                scatter(k,fxy(x,y),'red');
                xlabel('k');
                ylabel('f(xk,yk)');
                legend('f(x,y)','f(xk,yk)');
                hold off;

        end
    end
end

%% function to execute Newton method with gk constant
function [x,y] = newton_gk_const(x1,y1,gk,epsilon)
    x(1) = x1;
    y(1) = y1;
    k = 1;
    while norm(grad_fxy(x(k),y(k))) >= epsilon
        d(:,k) = -inv(grad2_fxy(x(k),y(k)))*grad_fxy(x(k),y(k));
        x(k+1) = x(k)+gk*d(1,k);
        y(k+1) = y(k)+gk*d(2,k);
        k = k + 1;       
    end
end

%% function to execute Newton method with gk that minimizes f(xk+gk*dk)
function [x,y] = newton_gk_minf(x1,y1,epsilon)
    x(1) = x1;
    y(1) = y1;
    k = 1;
    while norm(grad_fxy(x(k),y(k))) >= epsilon
        d(:,k) = -inv(grad2_fxy(x(k),y(k)))*grad_fxy(x(k),y(k));
        eqn = dfxy(x(k),y(k)) == 0;
        g_poss(:,1) = solve(eqn,'Real',true);
        g_poss(:,1) = vpa(g_poss(:,1));
        for i = 1:length(g_poss(:,1))
            F(i) = vpa(fxy(x(k)+g_poss(i,1)*d(1,k),y(k)+g_poss(i,1)*d(2,k)));
        end
        [~,I] = sort(F);
        g(1,k) = g_poss(I(1),1);
        x(k+1) = x(k)+g(1,k)*d(1,k);
        y(k+1) = y(k)+g(1,k)*d(2,k);
        k = k + 1;       
    end
end

%% function to execute Newton method with gk from armijo
function [x,y] = newton_gk_armijo(x1,y1,epsilon)
    x(1) = x1;
    y(1) = y1;
    k = 1;
    while norm(grad_fxy(x(k),y(k))) >= epsilon
        d(:,k) = -inv(grad2_fxy(x(k),y(k)))*grad_fxy(x(k),y(k));
        a = 0.01;
        b = 0.2;
        s = 1;
        m = 0;
        g(1,k) = s*b^m;
        x(k+1) = x(k)+g(1,k)*d(1,k);
        y(k+1) = y(k)+g(1,k)*d(2,k);
        l = -a*b^m*s*d(:,k).'*grad_fxy(x(k),y(k));
        while (fxy(x(k),y(k)) - fxy(x(k+1),y(k+1))) < l
            m = m + 1;
            g(k) = s*b^m;
            l = -a*b^m*s*d(:,k).'*grad_fxy(x(k),y(k));
            x(k+1) = x(k)+g(1,k)*d(1,k);
            y(k+1) = y(k)+g(1,k)*d(2,k);
        end
        x(k+1) = x(k)+g(1,k)*d(1,k);
        y(k+1) = y(k)+g(1,k)*d(2,k);
        k = k + 1;
        if(k > 50)
            fprintf('The newton method with gk from armijo is not working here (x1,y1) = (%.0f,%.0f)',x1,y1);
            break;
        end
    end
end
%% function to calculate the f(x,y)
function fxy = fxy(xk,yk)
    syms x y;
    fxy = x^5*exp(-x^2-y^2);
    fxy = subs(fxy,{x,y},{xk,yk});
end
%% function to calculate the gradient of f(x,y) at X,Y
function [grad_fXY] = grad_fxy(X,Y)
    syms x y;
    fxy = x^5*exp(-x^2-y^2);
    grad_fxy = gradient(fxy,[x,y]);
    grad_fXY = subs(grad_fxy,[x y],{X,Y});
end
%% function to calculate the 2nd grade gradient of f(x,y) at X,Y
function [grad2_fXY] = grad2_fxy(X,Y)
    syms x y;
    fxy = x^5*exp(-x^2-y^2);
    grad_fxy = gradient(fxy,[x,y]);
    grad2_fxy(:,1) = gradient(grad_fxy(1,1).',[x,y]);
    grad2_fxy(:,2) = gradient(grad_fxy(2,1).',[x,y]);
    grad2_fXY = subs(grad2_fxy,[x y],{X,Y});
end
%% function to calculate the df(xk - gk*grad(f(xk)))/dgk
function dfxy = dfxy(xk,yk)
    syms x y gk;
    fxy = x^5*exp(-x^2-y^2);
    d(:,1) = -grad_fxy(xk,yk);
    fxy = subs(fxy,{x,y},{xk+gk*d(1,1),yk+gk*d(2,1)});
    dfxy = diff(fxy,1,'gk');
end



