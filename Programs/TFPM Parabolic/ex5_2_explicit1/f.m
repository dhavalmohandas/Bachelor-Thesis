function [ a ] = f( x,t )

global Eps c

a = 2*exp(2*(x-1)/Eps)*cos(2*t) + exp((x-1))*cos(t)-Eps*exp((x-1))*sin(t)...
    +2*exp((x-1))*sin(t) + c*exp(2*(x-1)/Eps)*sin(2*t) + c*exp(x-1)*sin(t);

end

