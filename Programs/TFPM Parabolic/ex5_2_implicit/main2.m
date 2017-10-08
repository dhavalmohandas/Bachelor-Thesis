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
dt = 1e-4;
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
for i=1:N
    u0(i) = u0(i)-f(x(i),0)/c;
end
u = zeros(N,1);
A1 = zeros(N-2);
A2 = zeros(N-2);


for t=1:T
    for j=2:N-1
        a1 = (exp(-c*(j-1)*h*dt/2+b*c*dt^2*(t+1/4))-exp(-c*dt/2))*(1-exp(b*h/eps))/(exp(-c*t*dt*h-b*h/eps)-exp(-c*t*dt*h)-exp(-b*h/eps)-exp(b*h/eps+c*t*dt)+exp(b*h/eps)+exp(c*t*dt*h));
        a3 = (exp(-c*(j-1)*h*dt/2+b*c*dt^2*(t+1/4))-exp(-c*dt/2))*(exp(-b*h/eps)-1)/(exp(-c*t*dt*h-b*h/eps)-exp(-c*t*dt*h)-exp(-b*h/eps)-exp(b*h/eps+c*t*dt)+exp(b*h/eps)+exp(c*t*dt*h));
        a2 = exp(-c*dt/2)-a1-a3;
        b1 = (exp(c*dt*(j-1)*h/2-b*c*dt^2*(t+3/4))-exp(c*dt/2))*(1-exp(b*h/eps))/(exp(-c*(t+1)*dt*h-b*h/eps)-exp(-c*(t+1)*dt*h)-exp(-b*h/eps)-exp(b*h/eps+c*(t+1)*dt*h)+exp(b*h/eps)+exp(c*(t+1)*dt/h));
        b3 = (exp(c*dt*(j-1)*h/2-b*c*dt^2*(t+3/4))-exp(c*dt/2))*(exp(-b*h/eps)-1)/(exp(-c*(t+1)*dt*h-b*h/eps)-exp(-c*(t+1)*dt*h)-exp(-b*h/eps)-exp(b*h/eps+c*(t+1)*dt*h)+exp(b*h/eps)+exp(c*(t+1)*dt/h));
        b2 = exp(c*dt/2)-b1-b3;
        A1(j-1,j-1) = a2;
        A2(j-1,j-1) = b2;
        if j-1>1
            A1(j-1,j-2) = a1;
            A2(j-1,j-2) = b1;
        end
        if j<N-2
            A1(j-1,j) = a3;
            A2(j-1,j) = b3;
        end
        
    end
    %boundary conditions
    u(1) = exp(-2/eps)*sin(2*t*dt) + sin(t*dt)/exp(1) - f(x(1),t*dt)/c;
    u(N) = sin(2*t*dt) + sin(t*dt)-f(x(N),t*dt)/c;
    R = A1*u0(2:N-1);
    R(1) = R(1) + a1*u0(1)-b1*u(1);
    R(N-2) = R(N-2) + a3*u0(N)-b3*u(N);
    u(2:N-1) = A2\R;
    
    u0 = u;
end

v = zeros(N,1);

for j=1:N
    v(j) = u(j) + f(x(j),tf)/c;
end

%exact solution
for j=1:N
    uex(j,1) = exp(2*(x(j)-1)/eps)*sin(2*tf) + exp(x(j)-1)*sin(tf);
end

plot(x,v)
hold on
plot(x,uex)
legend('approximate solution','exact solution')
max(abs(v-uex))