clear all
close all

xi = 0;
xf = 2;
ti = 0;
tf = 1;
h  = 1/2^6;
dt = 1/2^4;
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
T = fix((tf-ti)/dt);
m = 4;
u0 = zeros(N,1);
u = zeros(N,1);
temp = zeros(N-2,1);

%initial condition

for j=1:N
    u0(j) = exp(-50*x(j)^2)*exp(1i*sin(x(j)));
end

k = cos(x);
A1 = zeros(N-2);
A2 = zeros(N-2);
a = 1;

for t=1:T
    for j=2:N-1
        a1(j) = sin(k(j)*(a*dt-(m+1)*h)/2)/sin(k(j)*(a*dt-(m-1)*h)/2);
        a2 = 1;
        A1(j-1,j-1) = 1;
        A2(j-1,j-1) = -a1(j);
        if j>2
            A1(j-1,j-2) = -a1(j);
            A2(j-1,j-2) = a2;
        end
        
    end
    
    for c=1:m
        temp(c) = u0(1);
    end
    
    temp(m+1:N-2) = u0(2:N-(m+1));
    u(1) = exp(-50*(x(1)-t*dt)^2)*exp(1i*sin(x(1)-t*dt));
    u(N) = exp(-50*(x(N)-t*dt)^2)*exp(1i*sin(x(N)-t*dt));
    R = A2*temp;
    R(1) = R(1) +a1(2)*u(1) + a2*u0(1);
    u(2:N-1) = A1\R;
    u0 = u;
end



usol = zeros(N,1);

for j=1:N
    usol(j) = exp(-50*(x(j)-tf)^2)*exp(1i*sin(x(j)-tf));
end

figure(1)
plot(x,real(usol))
hold on
plot(x,real(u))
title('real part')
legend('exact solution','approx solution')

figure(2)
plot(x,imag(usol))
hold on
plot(x,imag(u))
title('imaginary part')
legend('exact solution','approx solution')

max(abs(usol-u))