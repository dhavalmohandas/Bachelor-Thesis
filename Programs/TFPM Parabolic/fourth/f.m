function [ a ] = f( x,t )

global eps c

a = 2*exp(2*(x-1)/eps)*cos(2*t) + exp((x-1))*cos(t)-eps*exp((x-1))*sin(t)...
    +2*exp((x-1))*sin(t) + c*exp(2*(x-1)/eps)*sin(2*t) + c*exp(x-1)*sin(t);

end

