clear all
close all
clc

format long

xi = -5;
xf = 5;
ti = 0;
tf = 1;
dt = 0.05;
T = fix((tf-ti)/dt);
N = 20;
x = linspace(xi,xf,N);
u0 = zeros(1,N);
u = zeros(1,N);

%initial condition
for i=1:N
    u0(i) = exp(-5*x(i));           
end

for j=1:T
    
    for i=2:N-1
        alpha = local(i);
        u(i) = alpha(1)*u0(i-1)+alpha(2)*u0(i)*alpha(3)*u0(i+1);
    end
    u0 = u;
end

% exact solution
uexact = zeros(N,1);

for i=1:N
    uexact(i) = exp(-5*(x(i)-tf));
end

plot(x,u);
plot(x,uexact);


