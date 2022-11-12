%Lab 1 Task 1
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%initial interval
a1 = -1;
b1 = 3;

for case_of_fxi = 1:3 %loop that runs 3 times, each time for one of the three functions
    M1 = 1000; %length of the different lamdas or epsilons (depending on which is changeable)
    %first case: l = 0.01 (const) and e changeable
    l_c1 = 0.01 * ones(1,M1);
    e_c1 = linspace(0.00004,0.004,M1);
    [k_c1,a_c1,b_c1,calcs_c1] = bisector_method(case_of_fxi,l_c1,e_c1,a1,b1,M1);
    %second case: l changeable and e = 0.001 (const)
    l_c2 = linspace(0.005,0.05,M1); 
    e_c2 = 0.001 * ones(1,M1);
    [k_c2,a_c2,b_c2,calcs_c2] = bisector_method(case_of_fxi,l_c2,e_c2,a1,b1,M1);
    %third case: l = 0.01, e = 0.001
    l_c3 = 0.01;
    e_c3 = 0.001;
    M2 = 1;
    [k_c3,a_c3,b_c3,calcs_c3] = bisector_method(case_of_fxi,l_c3,e_c3,a1,b1,M2);

    %plots of the results for each case of study
    %first case: plot the calcualtions of fxi for l = 0.01 (const) and e changeable
    figure();
    plot(e_c1,calcs_c1,'k');
    title(['lamda = 0.01 and epsilon changeable in [0.00004,0.004] for f',num2str(case_of_fxi)]);
    xlabel('epsilon');
    ylabel('fxi calculations');
    %second case: plot the calcualtions of fxi for l changeable and e = 0.001 (const)
    figure();
    plot(l_c2,calcs_c2,'k');
    title(['lamda changeable in [0.005,0.05] and epsilon = 0.001 for f',num2str(case_of_fxi)]);
    xlabel('lamda');
    ylabel('fxi calculations');
    %third case: plor (k,ak) and (k,bk) for l = 0.01, e = 0.001
    figure()
    plot(a_c3,'b--o');
    title(['[ak,bk] as (k,ak) and (k,bk) for f',num2str(case_of_fxi)]);
    hold on;
    plot(b_c3,'r--o');
    xlabel('k');
    ylabel('[ak,bk]');
    legend('ak','bk');
    hold off;
end

%function to easily calculate the fxi depending on the case
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

%the implementation of the bisector method with landa and epsilon as vectors in order to be
%able to check the results with a constant value of either l or e and changeable value of the other
function [k,a,b,calcs] = bisector_method(case_of_fxi,l,e,a1,b1,M)
a(1)=a1; %I dont knpw the size of a so I will use dynamic memory allocation
b(1)=b1; %I dont knpw the size of b so I will use dynamic memory allocation
k = ones(1,M); % k is of size 1XM. M is the length of the different lamdas or epsilons 
calcs = zeros(1,M); %to store every fxi calculation
i = 1;
while (b(k(i)) - a(k(i))) >= l(i) %the loop terminates when I have checked that bk-ak < l for all possible lamdas (const or changeable)
    x1(k(i)) = (a(k(i))+b(k(i)))/2 - e(i);  %e can be const or changeable
    x2(k(i)) = (a(k(i))+b(k(i)))/2 + e(i);  %e can be const or changeable
    y1 = f(x1(k(i)),case_of_fxi);
    y2 = f(x2(k(i)),case_of_fxi);
    calcs(i) = calcs(i) + 2; %each time I calculate the value of fxi, I increase reps(i)
    if (y1 < y2) 
        a(k(i)+1) = a(k(i));
        b(k(i)+1) = x2(k(i));
    else 
        a(k(i)+1) = x1(k(i));
        b(k(i)+1) = b(k(i));
    end
    k(i) = k(i)+1;
    if ((b(k(i)) - a(k(i)) < l(i)) && (i < M)) %if bk-ak<l and I havent checked all lamda cases
        %when I have found the [ak,bk] with the x* for a spesific lamda, I increase i so I can continue with the next lamda
        i = i + 1;
    end    
end
end
    



    
