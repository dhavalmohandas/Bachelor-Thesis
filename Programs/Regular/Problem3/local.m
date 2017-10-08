function coeff = local(p)

xi = -5;
xf = 5;
N = 20;
a = 1;
dt = 0.05;

x = linspace(xi,xf,N);
h = (xf-xi)/(N-1);

k = -5*exp(-5*x(p));


alpha1 = sin(0.5*k*a*dt)*sin(0.5*k*(a*dt+h))/(sin(0.5*k*h)*sin(k*h));
alpha3 = sin(0.5*k*a)*sin(k*(0.5*a*dt-h))/(sin(0.5*k*h)*sin(k*h));
alpha2 = 1-alpha1-alpha3;

coeff = [alpha1,alpha2,alpha3];

end
