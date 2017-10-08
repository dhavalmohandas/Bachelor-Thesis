clear all
close all
clc
global eps c
xi = 0;
xf = 1;
ti = 0;
tf = 1;
b = 2;
c = 4;
dt = 2e-4;
eps = 0.01;
h = 0.1;
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
T = fix((tf-ti)/dt);
mu = dt/h^2;
A1 = zeros(N-2);
A2 = zeros(N-2);
%intial condition
u0 = zeros(N,1);
u = zeros(N,1);
v = zeros(N,1);

for i=1:N-2
    A1(i,i) = 10+12*mu;
    A2(i,i) = 10-12*mu;
    if i>1
        A1(i,i-1) = 1-6*mu;
        A2(i,i-1) = 1+6*mu;
    end
    if i<N-2
        A1(i,i+1) = 1-6*mu;
        A2(i,i+1) = 1+6*mu;
    end
end

for t=1:T
    u(1) = exp(2*(x(1)-1)/eps)*sin(2*t*dt) + exp(x(1)-1)*sin(t*dt);
    u(N) = exp(2*(x(N)-1)/eps)*sin(2*t*dt) + exp(x(N)-1)*sin(t*dt);
    R = A2*u0(2:N-1);
    R(1) = R(1) + (1+6*mu)*u0(1)-(1-6*mu)*u(1);
    R(N-2) = R(N-2) + (1+6*mu)*u0(N)-(1-6*mu)*u(N);
    u(2:N-1) = A1\R;
    u0 = u;
end

for i=1:N
    v(i) = u(i) + f(x(i),tf)/c;
end

%exact solution
for i=1:N
    uex(i,1) = exp(2*(x(i)-1)/eps)*sin(2*tf) + exp(x(i)-1)*sin(tf);
end

plot(x,v)
hold on
plot(x,uex)
legend('approx solution','exact solution')

max(abs(uex-v))