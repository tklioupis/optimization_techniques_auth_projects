%Lab 1 Task 2
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%initial interval
a1 = -1;
b1 = 3;

for case_of_fxi = 1:3 %loop that runs 3 times, each time for one of the three functions
    M1 = 1000; %length of the different lamdas or epsilons (depending on which is changeable)
    %first case: l changeable
    l_c1 = linspace(0.005,0.05,M1); 
    [k_c1,a_c1,b_c1,calcs_c1] = golden_sector_method(case_of_fxi,l_c1,a1,b1,M1);
    %second case: l = 0.01
    l_c2 = 0.01;
    M2 = 1;
    [k_c2,a_c2,b_c2,calcs_c2] = golden_sector_method(case_of_fxi,l_c2,a1,b1,M2);
    %third case: l = 0.02
    l_c3 = 0.02;
    M3 = 1;
    [k_c3,a_c3,b_c3,calcs_c3] = golden_sector_method(case_of_fxi,l_c3,a1,b1,M3);
    %fourth case: l = 0.03
    l_c4 = 0.04;
    M4 = 1;
    [k_c4,a_c4,b_c4,calcs_c4] = golden_sector_method(case_of_fxi,l_c4,a1,b1,M4);
    
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

%the implementation of the golden sector method with landa as vector in order to be
%able to check the results with a constant lamda value or changeable
function [k,a,b,calcs] = golden_sector_method(case_of_fxi,l,a1,b1,M)
g = 0.618;
a(1)=a1; %I dont knpw the size of a so I will use dynamic memory allocation
b(1)=b1; %I dont knpw the size of b so I will use dynamic memory allocation
k = ones(1,M); % k is of size 1XM. M is the length of the different lamdas or epsilons 
calcs = zeros(1,M); %to store every fxi calculation
i = 1;
while (b(k(i)) - a(k(i))) >= l(i) %the loop terminates when I have checked that bk-ak < l for all possible lamdas (const or changeable)
    if (k(i) == 1) %initialization
        x1(k(i)) = a(k(i)) + (1-g)*(b(k(i))-a(k(i)));  
        x2(k(i)) = a(k(i)) + g*(b(k(i))-a(k(i)));
        y1(k(i)) = f(x1(k(i)),case_of_fxi);
        y2(k(i)) = f(x2(k(i)),case_of_fxi);
        calcs(i) = calcs(i) + 2; %each time I calculate the value of fxi, I increase reps(i)
    end
    if (y1(k(i)) > y2(k(i))) 
        a(k(i)+1) = x1(k(i));
        b(k(i)+1) = b(k(i));

        x2(k(i)+1) = a(k(i)+1) + g*(b(k(i)+1) - a(k(i)+1));
        x1(k(i)+1) = x2(k(i));

        y1(k(i)+1) = y2(k(i));
        y2(k(i)+1) = f(x2(k(i)+1),case_of_fxi);
        calcs(i) = calcs(i) + 1;
    else 
        a(k(i)+1) = a(k(i));
        b(k(i)+1) = x2(k(i));

        x2(k(i)+1) = x1(k(i));
        x1(k(i)+1) = a(k(i)+1) + (1-g)*(b(k(i)+1) - a(k(i)+1));

        y1(k(i)+1) = f(x1(k(i)+1),case_of_fxi);
        y2(k(i)+1) = y1(k(i));
        calcs(i) = calcs(i) + 1;
    end
    k(i) = k(i)+1;
    if ((b(k(i)) - a(k(i)) < l(i)) && (i < M)) %if bk-ak<l and I havent checked all lamda cases
        %when I have found the [ak,bk] with the x* for a spesific lamda, I increase i so I can continue with the next lamda
        i = i + 1;
    end    
end
end


