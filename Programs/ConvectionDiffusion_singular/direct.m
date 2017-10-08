
close all
clear all
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
M = (N-2)^2;
A = zeros(M);
alpha = local(2,2)
v = zeros(N,N);
uapprox = zeros(N,N);

%boundary conditions
for i=1:N
    v(i,1) = (-f(x(i),y(1))/b)*exp(-(p*x(i)+q*y(1))/(2*epsilon^2));
    v(i,N) = (-f(x(i),y(N))/b)*exp(-(p*x(i)+q*y(N))/(2*epsilon^2));
    v(1,i) = (-f(x(1),y(i))/b)*exp(-(p*x(1)+q*y(i))/(2*epsilon^2));
    v(N,i) = (-f(x(N),y(i))/b)*exp(-(p*x(N)+q*y(i))/(2*epsilon^2));
end
v
% filling the matrix
for i=1:M
    A(i,i) = alpha(3);
    if i>4
        A(i,i-4) = alpha(4);
    end
    
    if i<M-3
        A(i,i+4) = alpha(5);
    end
    
    if i<M
        A(i,i+1) = alpha(2);
    end
    
    if i>1
        A(i,i-1) = alpha(1);
    end
end

i= N-2;

while i<M
    A(i+1,i) = 0;
    A(i,i+1) = 0;
    i = i+(N-2);    
end


% R.H.S of matrix equation
B = zeros(M,1);
k=1;
for i=2:(N-1)
    for j=2:(N-1)
        if i==2
            B(k) = B(k) -alpha(5)*v(i-1,j);
        end
        if i==N-2
            B(k) = B(k) -alpha(4)*v(i+1,j);            
        end
        if j==2
            B(k) = B(k) -alpha(2)*v(i,j-1);            
        end
        if j==N-2
            B(k) = B(k) -alpha(1)*v(i,j+1);            
        end
        k = k+1;
        
    end
   
end


temp = A\B;

k=1;

for i=2:N-1
    for j=2:N-1
        uapprox(i,j) = f(x(i),y(j))/b + temp(k)*exp( ( p*x(i)+q*y(j) )/( 2*epsilon^2 ) );
        k = k+1;
    end
end

usol = zeros(N,N);

for i=1:N
    for j=1:N
        usol(i,j) = y(j)*(1-y(j))*( exp((x(i)-1)/epsilon^2) + (x(i)-1)* exp(-1/epsilon^2) -x(i) );
    end
end

err = max(max(abs(usol-uapprox)))

figure(1);
surf(x,y,uapprox);
title('The approximate solution')

figure(2);
surf(x,y,usol);
title('The exact solution')
