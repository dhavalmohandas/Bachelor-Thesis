clear all
close all
clc
format long
xi = 0;
xf = 2;
epsilon = 0.2;
ti = 0;
tf = 1;
k = 1;
m = 2;
h = 1/2^7;
dt = 1/2^6;
T = fix((xf-xi)/dt);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
u = zeros(N,1);
u0 = zeros(N,1);
usol = zeros(N,1);

%initial condition

for j=1:N
    u0(j) = sign(x(j))*exp(-50*x(j)^2)*exp(1i*x(j)/epsilon);
end

%exact solution
for j=1:N
    usol(j) = sign(x(j)-tf)*exp(-50*(x(j)-tf)^2)*exp(1i*(x(j)-tf)/epsilon);
end


figure(1)
plot(x,real(usol))
xlabel('x')
ylabel('real(u)')
legend('exact solution')
title('real part')

figure(2)
plot(x,imag(usol))
xlabel('x')
ylabel('imag(u)')
legend('exact solution')
title('imaginary part')

