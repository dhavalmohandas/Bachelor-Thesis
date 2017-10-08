function [ a ] = f( x,t )
global Eps 
a = exp((x-1)/2)*(cos(t)-(Eps*sin(t))/4+sin(t)/2)+...
    exp((x-1)/Eps)*(2*cos(t)-(1/Eps)*sin(2*t) + (1/Eps)*sin(2*t));
%a = exp((x-1)/2)*(cos(t)-Eps*sin(t)/4-b*sin(t)/2)+...
%    exp((x-1)/eps)*(2*cos(2*t)-(1+b)*sin(2*t)/eps);
end

