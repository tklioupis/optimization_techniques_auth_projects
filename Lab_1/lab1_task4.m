%Lab 1 Task 4
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%% initial interval
a1 = -1;
b1 = 3;

%% run the method for each of the three functions and plot the results
for case_of_fxi = 1:3
    %first case: l changeable
    M1 = 100;   %length of the different lamdas for this case
                %ATTENTION: Because I use symbolic expressions do not use M > 100 because the script runs slow
    l_c1 = linspace(0.005,0.05,M1); 
    [k_c1,a_c1,b_c1,calcs_c1] = bisector_method_diff(case_of_fxi,l_c1,a1,b1,M1);
    %second case: l = 0.01
    l_c2 = 0.01;
    M2 = 1;
    [k_c2,a_c2,b_c2,calcs_c2] = bisector_method_diff(case_of_fxi,l_c2,a1,b1,M2);
    %third case: l = 0.02
    l_c3 = 0.02;
    M3 = 1;
    [k_c3,a_c3,b_c3,calcs_c3] = bisector_method_diff(case_of_fxi,l_c3,a1,b1,M3);
    %fourth case: l = 0.03
    l_c4 = 0.04;
    M4 = 1;
    [k_c4,a_c4,b_c4,calcs_c4] = bisector_method_diff(case_of_fxi,l_c4,a1,b1,M4);

    %plots
    %first case: l changeable , plot (l,calculations of fxi)
    figure();
    plot(l_c1,calcs_c1,'k');
    title(['lamda changeable in [0.005,0.05] for f',num2str(case_of_fxi)]);
    xlabel('lamda');
    ylabel('fxi calculations');
    %second,third and fourth case: plot (k,ak) and (k,bk) for l = 0.01, l = 0.02, l = 0.04
    figure()
    plot(a_c2,'b--o');
    title(['[ak,bk] as (k,ak) and (k,bk), lamda = [0.01,0.02,0.04], for f',num2str(case_of_fxi)]);
    hold on;
    plot(b_c2,'r--o');
    xline([k_c2 k_c3 k_c4],'-',{'l = 0.01','l = 0.02','l = 0.04'});
    xlabel('k');
    ylabel('[ak,bk]');
    legend('ak','bk');
    hold off;
end

%% the implementation of the bisector method with df(x)/dx with landa as a vector to check the results for multiple values of l
function [k,a,b,calcs] = bisector_method_diff(case_of_fxi,l,a1,b1,M)
    a(1) = a1; %I dont knpw the size of a so I will use dynamic memory allocation
    b(1) = b1; %I dont knpw the size of b so I will use dynamic memory allocation
    n(1) = ceil(log(l(1)/(b(1)-a(1)))/log(0.5));    %from theory n >= log0.5(l/(b1-a1))
                                                    %I use this expression logB(X) = logA(X) / logA(B)
    x(1) = (a(1)+b(1))/2;
    dfx = double(df(x(1),case_of_fxi));
    k = ones(1,M); % k is of size 1XM. M is the length of the different lamdas or epsilons 
    calcs = ones(1,M); %to store every fxi calculation, starts from 1
    i = 1;
    while (dfx ~= 0 && k(i) < n(i)) %the loop terminates when d(fx)/dx == 0 or n = k
        if dfx > 0 
            a(k(i)+1) = a(k(i));
            b(k(i)+1) = x(k(i));
        else 
            a(k(i)+1) = x(k(i));
            b(k(i)+1) = b(k(i));
        end
        k(i) = k(i) + 1;
        x(k(i)) = (a(k(i))+b(k(i)))/2;
        dfx = double(df(x(k(i)),case_of_fxi));
        calcs(i) = calcs(i) + 1;
        if ((dfx == 0 || k(i) >= n(i)) && (i < M))
            %when I have found the [ak,bk] with the x* for a spesific lamda, I increase i so I can continue with the next lamda
            i = i + 1;
            n(i) = ceil(log(l(i)/(b(1)-a(1)))/log(0.5));
        end    
    end
end

%% function to easily calculate the dfxi depending on the case
function [dfxi] = df(x, i)
    syms z;
    switch i
        case 1
            fxi = (z-2)^2 + z*log(z+3);
            dfxi = diff(fxi);
        case 2  
            fxi = 5^z + (2-cos(z))^2;
            dfxi = diff(fxi);
        case 3
            fxi = exp(z)*(z^3-1) + (z-1)*sin(z);
            dfxi = diff(fxi);
    end
    dfxi = subs(dfxi, z, x);
end