clear all
close all
clc
format long
xi = 0;
xf = 2;
h = 1/(2^6);
ti = 0;
tf = 1;
ti = 0;
tf = 1;
dt = 1/(2^7);

for p=1:4
T = fix((tf-ti)/dt);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
u = zeros(1,N);
%initial condition
for j=1:N
    u0(j) = exp(-50*x(j)^2)*exp(i*sin(x(j)));
end

k = cos(x);

for m=1:T
    
    for j=2:N-1
        a1 = sin(k(j)*dt/2)*sin(k(j)*(dt+h)/2)/(sin(k(j)*h/2)*sin(k(j)*h));
        a3 = sin(k(j)*dt/2)*sin(k(j)*(dt-h)/2)/(sin(k(j)*h/2)*sin(k(j)*h));
        a2 = 1-a1-a3;
        u(j) = a1*u0(j-1)+a2*u0(j)+a3*u0(j+1);
    end
    u(1) = exp(-50*(x(1)-m*dt)^2)*exp(i*sin(x(1)-m*dt));
    u(N) = exp(-50*(x(N)-m*dt)^2)*exp(i*sin(x(N)-m*dt));
    u0 = u;
end

uex = zeros(1,N); %exact solution

for j=1:N
    uex(j) = exp(-50*(x(j)-tf)^2)*exp(i*sin(x(j)-tf));
end
err(p) = max(abs(u-uex));
h = h/2;
dt = dt/2;
end

for i=1:p-1
    order(i) = log(err(i)/err(i+1))/log(2);
end


figure(1)
plot(x,real(u),'*')
hold on
plot(x,real(uex))
xlabel('x')
ylabel('real(u)')
title('real part')

figure(2)
plot(x,imag(u),'*')
hold on
plot(x,imag(uex))
xlabel('x')
ylabel('imag(u)')
title('imaginary part')

err
order