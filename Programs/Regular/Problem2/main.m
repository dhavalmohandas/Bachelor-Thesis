close all; 
clear all;
clc

N = 5;
tol = 1e-3;
xi = 0;
xf = 2;
yi = 0;
yf = 1;
x = linspace(xi,xf,2*N);
y = linspace(yi,yf,N);

u2 = zeros(2*N,N);

for i = 1:2*N
    u2(i,1) = 1; 
    u2(i,N) = exp(x(i));
end

for i=1:N
    u2(1,i) = 1; 
    u2(2*N,i) = exp(2*y(i));
end


err = 1000;
k = 0;
u = u2;

while err > tol
    for i = 2:2*N-1
        for j = 2:N-1
            alpha = local(i,j);
            u(i,j) = (x(i)^2+y(j)^2)*exp(x(i)*y(j))-( alpha(1)*u(i-1,j) + alpha(2)*u(i+1,j) + alpha(4)*u(i,j-1) + alpha(5)*u(i,j+1) )/alpha(3);
        end
       
    end
    err = max(max(abs(u-u2)))
    u2 = u;
    k = k+1
end

figure(1);
mesh(x,y,u2);
title('The approximate solution')

for i=1:2*N
    for j=1:N
        usol(i,j)=exp(x(i)*y(j));
    end
end
%{
figure(2);
mesh(x,y,usol);
title('The exact solution')
%}
error = max(max(abs(u-usol)))

