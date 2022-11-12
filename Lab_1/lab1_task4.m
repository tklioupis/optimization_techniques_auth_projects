%Lab 1 Task 4
%Theodoros Lioupis AEM 9733

clear all;
close all;
clc;

%initial interval
a1 = -1;
b1 = 3;

for case_of_fxi = 1:3

end

function [k,a,b,calcs] = bisector_method_diff(case_of_fxi,l,a1,b1,M)

end

result = double(df(2,4)); %example

%function to easily calculate the dfxi depending on the case
function [dfxi] = df(x, i)
    syms z;
    switch i
        case 1
            fxi = (z-2)^2 + z*log(z+3);
            dfxi = diff(fxi);
        case 2  
            fxi = 5^zeros + (2-cos(z))^2;
            dfxi = diff(fxi);
        case 3
            fxi = exp(z)*(z^3-1) + (z-1)*sin(z);
            dfxi = diff(fxi);
        case 4 %test
            fxi = 5*z^2;
            dfxi = diff(fxi);
        dfxi = subs(dfxi, z, x);
    end       
end