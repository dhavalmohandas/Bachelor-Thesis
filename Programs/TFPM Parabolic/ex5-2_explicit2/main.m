clear all
close all
clc
format long
global Eps c
xi = 0;
xf = 1;
ti = 0;
tf = 1;
b = 2;
c = 0;
h = 1/16;
dt = 2e-5;
Eps = h*0.01;

%a1 = b*dt/(h*(1-exp(-b*h/Eps)));
%a3 = exp(-b*h/Eps)*a1;
a1 = (b*dt/h)*(exp(b*h/Eps)-1)/(exp(b*h/Eps)+exp(-b*h/Eps)-2);
a3 = (b*dt/h)*(1-exp(-b*h/Eps))/(exp(b*h/Eps)+exp(-b*h/Eps)-2);
a2 = 1-a1-a3;

T = fix((tf-ti)/dt);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);

%initial condition
v0 = zeros(N,1);
%exact solution
uex = zeros(N,1);
v = zeros(N,1);

for t=1:T
    for i=2:N-1
        v(i) = a1*v0(i-1)+a2*v0(i)+a3*v0(i+1)+dt*f(x(i),t*dt);
    end
    %boundary conditions
    v(1) = exp(2*(x(1)-1)/Eps)*sin(2*t*dt) + exp(x(1)-1)*sin(t*dt);
    v(N) = exp(2*(x(N)-1)/Eps)*sin(2*t*dt) + exp(x(N)-1)*sin(t*dt);
    v0 = v;
end

u = zeros(N,1);
for i=1:N
    u(i) = v(i) + tf*f(x(i),tf);
end

%exact solution
for i=1:N
    uex(i,1) = exp(2*(x(i)-1)/Eps)*sin(2*tf) + exp(x(i)-1)*sin(tf);
end

plot(x,u)
hold on
plot(x,uex)
ylabel('u')
xlabel('x')
legend('approximate solution','exact solution')
max(abs(u-uex))
