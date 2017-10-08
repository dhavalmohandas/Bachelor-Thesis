clear all
close all
clc
format long
xi = 0;
xf = 2;
ti = 0;
tf = 1;
h = 1/(2^6);
dt = 1/(2^11);
T = fix((tf-ti)/dt);

for p=1:4

N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);



u = zeros(1,N);
u0 = zeros(1,N);

%initial condition
for j=1:N
    u0(j) = exp(-50*x(j)^2)*exp(1i*sin(x(j)));
end

k = cos(x);
a = 1;

for m=1:T
    u(1) = exp(-50*(x(1)-a*m*dt)^2)*exp(1i*sin(x(1)-a*m*dt));
    u(N) = exp(-50*(x(N)-a*m*dt)^2)*exp(1i*sin(x(N)-a*m*dt));
    for j=2:N-1
        a1 = sin(k(j)*a*dt/2)*sin(k(j)*(a*dt-h)/2)/(sin(k(j)*h/2)*sin(k(j)*h));
        a2 = -sin(k(j)*a*dt/2)*sin(k(j)*(a*dt-2*h)/2)/(sin(k(j)*h/2))^2;
        a3 = 1-a1-a2;
        if j==2
            u(j) = a1*u0(j-1)+a2*u0(j-1)+a3*u0(j);
        else
            u(j) = a1*u0(j-2)+a2*u0(j-1)+a3*u0(j);
        end
    end

    u0 = u;
end


uex = zeros(1,N); %exact solution

for j=1:N
    uex(j) = exp(-50*(x(j)-a*tf)^2)*exp(1i*sin(x(j)-a*tf));
end

h = h/2;
err(p) = max(abs(u-uex));

end

for i=1:p-1
    order(i) = log(err(i)/err(i+1))/log(2);
end

figure(1)
plot(x,real(u),'*')
hold on
plot(x,real(uex))
title('real part')
xlabel('x')
ylabel('real(u)')
legend('approx solution','exact solution')


figure(2)
plot(x,imag(u),'*')
hold on
plot(x,imag(uex))
title('imaginary part')
xlabel('x')
ylabel('imag(u)')
legend('approx solution','exact solution')

order