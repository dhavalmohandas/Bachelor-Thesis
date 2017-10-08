function [ a ] = f( x,t )

if x<=0.5;
    a = -2*x*t;
else
    a = 2*(x-1)*t;
end

end

