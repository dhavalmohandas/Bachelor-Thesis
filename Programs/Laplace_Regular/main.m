close all; 
clear all;
clc

N = 20;
tol = 1e-6;
xi = 0;
xf = 1;
yi = 0;
yf = 1;
x = linspace(xi,xf,N);
y = linspace(yi,yf,N);

u2 = zeros(N,N);

for i = 1:N
    u2(1,i) = -y(i)^2;
    u2(N,i) = 1-y(i)^2; 
    u2(i,1) = x(i)^2; 
    u2(i,N) = x(i)^2-1; 
end

err = 1000;
k = 0;
u = u2;
alpha = local(2,2);

while err > tol
    for i = 2:N-1
        for j = 2:N-1
            
            u(i,j) = -( alpha(1)*u(i-1,j) + alpha(2)*u(i+1,j) + alpha(4)*u(i,j-1) + alpha(5)*u(i,j+1) )/alpha(3);
        end
       
    end
    err = max(max(abs(u-u2)));
    u2 = u;
    
end

usol = zeros(N,N);

for i=1:N
    for j=1:N
        usol(i,j) = x(i)^2-y(j)^2;
    end
end

figure(1);
mesh(x,y,u2);
title('The approximate solution')

figure(2)
mesh(x,y,usol);
title('exact solution')

max(max(abs(usol-u2)))