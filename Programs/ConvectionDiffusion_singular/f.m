function a = f(x,y)
global epsilon

a= (2*epsilon^2 + y*(1-y))*(exp((x-1)/epsilon^2) + (x-1)*exp(-1/epsilon^2) + -x) + y*(1-y)*(exp(-1/epsilon^2)-1);

end