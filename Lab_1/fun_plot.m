%Lab 1
%Theodoros Lioupis AEM 9733

%plot the functions to see them and compare the results for the interval where the minimum is

for case_of_fxi = 1:3
    x = linspace(-1,3,500);
    y = f(x,case_of_fxi);
    figure()
    plot(x,y);
    title(['f',num2str(case_of_fxi)]);
end

%function to easily calculate the fxi depending on the case
function [fxi] = f(x, i)
    switch i
        case 1
            fxi = (x-2).^2 + x.*log(x+3);
        case 2  
            fxi = 5.^x + (2-cos(x)).^2;
        case 3
            fxi = exp(x).*(x.^3-1) + (x-1).*sin(x);
    end       
end