clear all
close all
format long

xi = 0;
xf = 1;
yi = 0;
yf = 1;
tol = 1e-3;
h = 1/64;
epsilon = 1 ;
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
y = linspace(xi,xf,N);


%boundary conditions
for i=1:N
    u(i,1) = 0;
    u(i,N) = 0;
    u(1,i) = 0;
    u(N,i) = 0;
end

u1 = u;
err = 100;

while err > tol
    for i=2:N-1
        for j=2:N-1
            b = x(i)+y(j)+1;
            p = x(i);
            q = y(j);
            mu = sqrt(b + (p^2+q^2)/(4*epsilon^2))/epsilon;
            T1 = (exp(-p*h/(2*epsilon^2))*u(i+1,j)+exp(-q*h/(2*epsilon^2))*u(i,j+1)+exp(p*h/(2*epsilon^2))*u(i-1,j)...
                +exp(q*h/(2*epsilon^2))*u(i,j-1))/(4*cosh(mu*h/2)^2);
            T2 = (f(x(i),y(j))/b)*(1-(exp(-p*h/(2*epsilon^2))+exp(-q*h/(2*epsilon^2))+exp(p*h/(2*epsilon^2))...
                +exp(q*h/(2*epsilon^2)))/(4*cosh(mu*h/2)^2));
            u(i,j) = T1+T2;
        end
        err = max(max(abs(u-u1)));
        u1 = u;
    end
end


surf(x,y,u1);