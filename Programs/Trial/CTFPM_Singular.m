clear all
close all
clc
format long
xi = 0;
xf = 1;
h = 1/(2^9);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
ti = 0;
tf = 0.25;
dt = 1/(2^10);
T = fix((tf-ti)/dt);
u = zeros(1,N);
epsilon = 0.01;

%initial condition
for j=1:N
    u0(j) = exp(-200*x(j)^2)*exp(i*sin(x(j))/epsilon);
end


a = 0.5+x;
%k = cos(x);

for m=1:T
    
    for j=2:N-1
        y = (0.5 + x(j))*exp(-m*dt)-0.5;
        k(j) = cos(y)*exp(-m*dt)/epsilon; 
        a1 = sin(k(j)*a(j)*dt/2)*sin(k(j)*(a(j)*dt+h)/2)/(sin(k(j)*h/2)*sin(k(j)*h));
        a3 = sin(k(j)*a(j)*dt/2)*sin(k(j)*(a(j)*dt-h)/2)/(sin(k(j)*h/2)*sin(k(j)*h));
        a2 = 1-a1-a3;
        u(j) = a1*u0(j-1)+a2*u0(j)+a3*u0(j+1);
    end
    u(1) = exp(-200*(x(1)-m*dt)^2)*exp(i*sin(x(1)-m*dt)/epsilon);
    u(N) = exp(-200*(x(N)-m*dt)^2)*exp(i*sin(x(N)-m*dt)/epsilon);
    u0 = u;
end


uex = zeros(1,N); %exact solution

for j=1:N
    uex(j) = exp(-200*(x(j)-tf)^2)*exp(i*sin(x(j)-tf)/epsilon);
end

figure(1)
plot(x,real(u))
hold on
plot(x,real(uex))
title('real part')

figure(2)
plot(x,imag(u))
hold on
plot(x,imag(uex))
title('imaginary part')

max(abs(u-uex))