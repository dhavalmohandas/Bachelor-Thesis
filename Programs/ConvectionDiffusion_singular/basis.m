function [a] = basis(i,x,y)

global mu

switch(i)
    case 1
        a = exp(-mu*x);
        return;
    case 2
        a = exp(mu*x);
        return;
    case 3
        a = exp(-mu*y);
        return;
    case 4
        a = exp(mu*y);
        return;
end

end
