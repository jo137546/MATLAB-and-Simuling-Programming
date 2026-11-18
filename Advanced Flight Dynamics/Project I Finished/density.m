%%
%Readme
%Written by: Elliott Eggers
%1/14/2023

% please check to see if this densicty model is correct or not}.

function rho = density(h)

if h<11000
    rho = 4.1748*10^-11*(288.14-0.00649*h)^4.256;
elseif h<25000
    rho = 0.3644*exp(1.73-0.000157*h);
else
    T = -131.21+0.00299*h;
    P = 2.488*((T+273.15)/216.6)^-11.388;
    rho = P/(0.2896*(T+273.15));
end