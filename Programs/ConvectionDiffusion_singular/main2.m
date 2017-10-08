close all; 
clear all;
%clc

N = 33;
tol = 1e-6;
xi = 0;
xf = 1;
yi = 0;
yf = 1;
x = linspace(xi,xf,N);
y = linspace(yi,yf,N);
h = (xf-xi)/(N-1);
global p q b epsilon
epsilon = 0.1;
p = 1;
q = 0;
b = 1;

mu = sqrt(b+(p^2+q^2)/(4*epsilon^2))/epsilon;
v = zeros(N,N);

% boundary condition
for i=1:N
    v(i,1) = 0;
    v(i,N) = 0;
    v(1,i) = 0;
    v(N,1) = 0;
end

v1 = v;

m = p*h/(2*epsilon^2);
err = 1000;

while err > tol
    for i=2:N-1
        for j=2:N-1
             v(i,j) = (f(x(i),y(j))/b)*(1-(1+cosh(m))/(2*cosh(mu*h*0.5)^2)) + ...
                 (exp(-m)*v(i+1,j) + v(i,j+1) + exp(m)*v(i-1,j) + v(i,j-1))/(4*cosh(mu*h*0.5)^2);
             %v(i,j) = ex(-x(i))*(v(i+1,j)*ex(x(i+1)) + v(i-1,j)*ex(x(i-1))+ex(x(i))*(v(i,j+1)+v(i,j-1) ))/(4*cosh(0.5*mu*h)^2)+...
                 %f(x(i),y(j))*((1-ex(-x(i))*( ex(x(i+1)) + ex(x(i-1)) + 2*ex(x(i)) )/(4*cosh(mu*0.5*h)^2)));
                 
        end
        err = max(max(abs(v-v1)));
        v1 = v;
    end
end



usol = zeros(N,N);

for i=1:N
    for j=1:N
        usol(i,j) = y(j)*(1-y(j))*( exp((x(i)-1)/epsilon^2) + (x(i)-1)* exp(-1/epsilon^2) -x(i) );
    end
end

figure(1);
surf(x,y,v1);
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
title('The approximate solution')

figure(2);
surf(x,y,usol);
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
title('The exact solution')
%max(max(abs(usol)))
max(max(abs(usol-v1)))