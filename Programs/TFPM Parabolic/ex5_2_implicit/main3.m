clear all
close all
format long
%clc
global eps c
xi = 0;
xf = 1;
ti = 0;
tf = 1;
b = 2;
c = 2;
eps = 0.1;  
h = 1/32;
dt = 2e-4;
%{
a1 = (exp(b*h/eps-c*dt/2)*b*dt)/(2*h*(exp(b*h/eps)-1));
a3 = (exp(-c*dt/2)*b*dt)/(2*h*(exp(b*h/eps)-1));
a2 = exp(-c*dt/2)-a1-a3;
b1 = -(exp(b*h/eps+c*dt/2)*b*dt)/(2*h*(exp(b*h/eps)-1));
b3 = -(exp(c*dt/2)*b*dt)/(2*h*(exp(b*h/eps)-1));
b2 = exp(c*dt/2)-b1-b3;
%}

a1 = exp(-c*dt/2)*(b*dt/(2*h))*(exp(b*h/eps)-1)/(exp(b*h/eps)+exp(-b*h/eps)-2);
a3 = exp(-c*dt/2)*(b*dt/(2*h))*(1-exp(-b*h/eps))/(exp(b*h/eps)+exp(-b*h/eps)-2);
a2 = exp(-c*dt)-a1-a3;

b1 = exp(c*dt/2)*(b*dt/(2*h))*(1-exp(b*h/eps))/(exp(-b*h/eps)+exp(b*h/eps)-2);
b3 = exp(c*dt/2)*(b*dt/(2*h))*(exp(-b*h/eps)-1)/(exp(-b*h/eps)+exp(b*h/eps)-2);
b2 = exp(c*dt/2)-b1-b3;

T = fix((tf-ti)/dt);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);

%initial condition
u0 = zeros(N,1);

u = zeros(N,1);
A1 = zeros(N-2);
A2 = zeros(N-2);

for i=1:N-2
    A1(i,i) = a2;
    A2(i,i) = b2;
    if i>1
        A1(i,i-1) = a1;
        A2(i,i-1) = b1;
    end
    if i<N-2
        A1(i,i+1) = a3;
        A2(i,i+1) = b3;
    end
end

for t=1:T
    %boundary conditions
    for j=2:N-1
        F(j-1,1) = ( b1*f(x(j-1),t*dt)-a1*f(x(j-1),(t-1)*dt)+b2*f(x(j),t*dt)-a2*f(x(j),(t-1)*dt)+b3*f(x(j+1),t*dt)-a3*f(x(j+1),(t-1)*dt) )/c;
    end
    u(1) = exp(-2/eps)*sin(2*t*dt) + sin(t*dt)/exp(1);
    u(N) = sin(2*t*dt) + sin(t*dt);
    R = A1*u0(2:N-1)+F;
    R(1) = R(1) + a1*u0(1)-b1*u(1);
    R(N-2) = R(N-2) + a3*u0(N)-b3*u(N);
    u(2:N-1) = A2\R;
    
    u0 = u;
end


%exact solution
for i=1:N
    uex(i,1) = exp(2*(x(i)-1)/eps)*sin(2*tf) + exp(x(i)-1)*sin(tf);
end

plot(x,u)
hold on
plot(x,uex)
legend('approximate solution','exact solution')
max(abs(u-uex))