close all; 
clear all;
clc

N = 12;
tol = 1e-6;
xi = 0;
xf = 1;
yi = 0;
yf = 1;
x = linspace(xi,xf,N);
y = linspace(yi,yf,N);
epsilon = 1;
p = 1;
q = 0;
b = 1;

u2 = zeros(N,N);

for i=1:N
    u2(i,1) = (-f(x(i),y(1))/b)*exp(-(p*x(i)+q*y(1))/(2*epsilon^2));
    u2(i,N) = (-f(x(i),y(N))/b)*exp(-(p*x(i)+q*y(N))/(2*epsilon^2));
    u2(1,i) = (-f(x(1),y(i))/b)*exp(-(p*x(1)+q*y(i))/(2*epsilon^2));
    u2(N,i) = (-f(x(N),y(i))/b)*exp(-(p*x(N)+q*y(i))/(2*epsilon^2));
end

uapprox = zeros(N,N);
err = 1000;
u = u2;
alpha = local(2,2);



while err > tol
    for i = 2:N-1
        for j = 2:N-1            
            u(i,j) = -( alpha(1)*u(i+1,j) + alpha(3)*u(i-1,j) + alpha(4)*u(i,j-1) + alpha(2)*u(i,j+1) )/alpha(3);
        end
       
    end
    err = max(max(abs(u-u2)));
    u2 = u;
    
end

for i=2:N-1
    for j=2:N-1
        uapprox(i,j) = f(x(i),y(j))/b + u2(i,j)*exp( ( p*x(i)+q*y(j) )/( 2*epsilon^2 ) );
    end
end


usol = zeros(N,N);

for i=1:N
    for j=1:N
        usol(i,j) = y(j)*(1-y(j))*( exp((x(i)-1)/epsilon^2) + (x(i)-1)* exp(-1/epsilon^2) -x(i) );
    end
end

figure(1);
surf(x,y,uapprox);
title('The approximate solution')

figure(2);
surf(x,y,usol);
title('The exact solution')

figure(3)
surf(x,y,u2);
title('trial')
%}
err = max(max(abs(uapprox-usol)))
err2 = max(max(abs(u2-usol)))