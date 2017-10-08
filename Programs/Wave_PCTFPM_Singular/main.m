clear all
close all
clc
format long
xi = 0;
xf = 2;
h = 1/(2^6);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
ti = 0;
tf = 1;
dt = 1/(2^4);
T = fix((tf-ti)/dt);
u = zeros(1,N);
epsilon = 0.001;

%initial condition
for j=1:N
    u0(j) = exp(-50*x(j)^2)*exp(i*sin(x(j))/epsilon);
end

k = 1/epsilon;
m = 4;

for p=1:T
    
    a1 = sin(k*(dt-(m+1)*h)/2)/sin(k*(dt-(m-1)*h)/2);
    b1 = 1;
    b0 = -a1;
    
    
    u(1) = exp(-50*(x(1)-p*dt)^2)*exp(i*sin(x(1)-p*dt)/epsilon);  %boundary conditions
    u(N) = exp(-50*(x(N)-p*dt)^2)*exp(i*sin(x(N)-p*dt)/epsilon);
    u0 = u;
end


uex = zeros(1,N); %exact solution

for j=1:N
    uex(j) = exp(-50*(x(j)-1)^2)*exp(i*sin(x(j)-1)/epsilon);
end

figure(1)
plot(x,real(u),'*')
hold on
plot(x,real(uex))
title('real part')

figure(2)
plot(x,imag(u),'*')
hold on
plot(x,imag(uex))
title('imaginary part')

max(abs(u-uex))