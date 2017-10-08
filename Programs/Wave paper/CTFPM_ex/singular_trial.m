
clear all
close all
clc
format long
xi = 0;
xf = 2;
h = 1/(2^8);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
ti = 0;
tf = 1;
ti = 0;
tf = 1;
dt = 1/(2^10);
T = fix((tf-ti)/dt);
u = zeros(1,N);
epsilon = 0.01;

%initial condition
for j=1:N
    u0(j) = exp(-50*x(j)^2)*exp(i*sin(x(j))/epsilon);
end

k = cos(x);


for m=1:T
    
    for j=2:N-1
        a1 = sin(k(j)*dt/2)*sin(k(j)*(dt+h)/2)/(sin(k(j)*h/2)*sin(k(j)*h));
        a3 = sin(k(j)*dt/2)*sin(k(j)*(dt-h)/2)/(sin(k(j)*h/2)*sin(k(j)*h));
        a2 = 1-a1-a3;
        u(j) = a1*u0(j-1)+a2*u0(j)+a3*u0(j+1);
    end
    u(1) = exp(-50*(x(1)-m*dt)^2)*exp(i*sin(x(1)-m*dt)/epsilon);
    u(N) = exp(-50*(x(N)-m*dt)^2)*exp(i*sin(x(N)-m*dt)/epsilon);
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