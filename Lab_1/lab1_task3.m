%Lab 1 Task 3
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%% initial interval
a1 = -1;
b1 = 3;

for case_of_fxi = 1:3 %loop that runs 3 times, each time for one of the three functions
    %first case: l changeable
    M1 = 100; %length of the different lamdas for this case
    l_c1 = linspace(0.005,0.05,M1); 
    e_c1 = 0.001;
    [k_c1,a_c1,b_c1,calcs_c1] = fibonacci_method(case_of_fxi,l_c1,e_c1,a1,b1,M1);
    %second case: l = 0.01
    l_c2 = 0.01;
    e_c2 = 0.001;
    M2 = 1;
    [k_c2,a_c2,b_c2,calcs_c2] = fibonacci_method(case_of_fxi,l_c2,e_c2,a1,b1,M2);
    %third case: l = 0.02
    l_c3 = 0.02;
    e_c3 = 0.001;
    M3 = 1;
    [k_c3,a_c3,b_c3,calcs_c3] = fibonacci_method(case_of_fxi,l_c3,e_c3,a1,b1,M3);
    %fourth case: l = 0.03
    l_c4 = 0.04;
    e_c4 = 0.001;
    M4 = 1;
    [k_c4,a_c4,b_c4,calcs_c4] = fibonacci_method(case_of_fxi,l_c4,e_c4,a1,b1,M4);
    
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


%% the implementation of the fibonacci method with landa as vector in order to be able to check the results with a constant lamda value or changeable
function [k,a,b,calcs] = fibonacci_method(case_of_fxi,l,e,a1,b1,M)
    a(1) = a1; %I dont knpw the size of a so I will use dynamic memory allocation
    b(1) = b1; %I dont knpw the size of b so I will use dynamic memory allocation
    k = ones(1,M); % k is of size 1XM. M is the length of the different lamdas or epsilons 
    calcs = zeros(1,M); %to store every fxi calculation
    %calculating fibonacci numbers
    fibb = ones(1,1000);
    %fibb(0) = 1;
    %fibb(1) = 1;
    for i = 3:1000
        fibb(i) = fibb(i-1) + fibb(i-2);
    end
    %finding the calculations of fxi for every lamda
    for i = 1:M % i -> case of lamda
        for j = 1:1000 % j -> which fibonacci number I am checking
            if fibb(j) > ((b(1) - a(1))/l(i))
                calcs(i) = j;
                break;
            end
        end
    end
    i = 1;
    while (k(i) < (calcs(i) - 2)) %the loop terminates when k >= (n - 2) for all possible lamdas (const or changeable)       
        if (k(i) == 1) %initialization
            x1(1) = a(1) + fibb(calcs(i)-2)/fibb(calcs(i))*(b(1)-a(1));
            x2(1) = a(1) + fibb(calcs(i)-1)/fibb(calcs(i))*(b(1)-a(1));
            y1(1) = f(x1(1),case_of_fxi);
            y2(1) = f(x2(1),case_of_fxi);
        end
        if (y1(k(i)) > y2(k(i))) %step 2
            a(k(i)+1) = x1(k(i));
            b(k(i)+1) = b(k(i));
            x1(k(i)+1) = x2(k(i));
            x2(k(i)+1) = a(k(i)+1) + fibb(calcs(i)-k(i)-1)/fibb(calcs(i)-k(i))*(b(k(i)+1) - a(k(i)+1));
            if(k == (calcs(i) - 2)) %step 5
                x1(calcs(i)) = x1(calcs(i)-1);
                x2(calcs(i)) = x1(calcs(i)-1) + e;
                y1(calcs(i)) = f(x1(calcs(i)),case_of_fxi);
                y2(calcs(i)) = f(x2(calcs(i)),case_of_fxi);
                if (y1(calcs(i)) > y2(calcs(i)))
                    a(calcs(i)) = x1(calcs(i));
                    b(calcs(i)) = b(calcs(i)-1);
                else
                    a(calcs(i)) = a(calcs(i)-1);
                    b(calcs(i)) = x2(calcs(i));
                end
                %algorithm must stop, moving to next lamda
                if(i<M)
                    i = i + 1;
                end
            else %step 4
                y1(k(i)+1) = f(x1(k(i)+1),case_of_fxi);
                y2(k(i)+1) = f(x2(k(i)+1),case_of_fxi);
                k(i) = k(i) + 1;
            end
        else %step 3
            a(k(i)+1) = a(k(i));
            b(k(i)+1) = x2(k(i));
            x2(k(i)+1) = x1(k(i));
            x1(k(i)+1) = a(k(i)+1) + fibb(calcs(i)-k(i)-2)/fibb(calcs(i)-k(i))*(b(k(i)+1) - a(k(i)+1));
            if(k == (calcs(i) - 2)) %step 5
                x1(calcs(i)) = x1(calcs(i)-1);
                x2(calcs(i)) = x1(calcs(i)-1) + e;
                y1(calcs(i)) = f(x1(calcs(i)),case_of_fxi);
                y2(calcs(i)) = f(x2(calcs(i)),case_of_fxi);
                if (y1(calcs(i)) > y2(calcs(i)))
                    a(calcs(i)) = x1(calcs(i));
                    b(calcs(i)) = b(calcs(i)-1);
                else
                    a(calcs(i)) = a(calcs(i)-1);
                    b(calcs(i)) = x2(calcs(i));
                end
                %algorithm must stop, moving to next lamda
                if(i<M)
                    i = i + 1;
                end
            else %step 4
                y1(k(i)+1) = f(x1(k(i)+1),case_of_fxi);
                y2(k(i)+1) = f(x2(k(i)+1),case_of_fxi);
                k(i) = k(i) + 1;
            end
        end    
    end
end

%% function to easily calculate the fxi depending on the case
function [fxi] = f(x, i)
    switch i
        case 1
            fxi = (x-2)^2 + x*log(x+3);
        case 2  
            fxi = 5^x + (2-cos(x))^2;
        case 3
            fxi = exp(x)*(x^3-1) + (x-1)*sin(x);
    end       
end