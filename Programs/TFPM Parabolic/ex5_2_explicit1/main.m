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
c = 1;
dt = 2e-5;
h = 1/16;
%Eps = c*h*1e-3;
%Eps = 1e-4;
Eps = 0.1;
L1 = b/(2*Eps)+sqrt(b^2/(4*Eps^2)+c/Eps);
L2 = b/(2*Eps)-sqrt(b^2/(4*Eps^2)+c/Eps);


for p=1:4    


%a1 = exp(-L1*h)*(1-exp(-c*dt))/((1-exp(-L1*h))*(1-exp(L2*h)));
a1 = (1-exp(-c*dt))*(exp(L2*h)-exp(L1*h))/(exp(-L1*h)*exp(L2*h)-exp(-L1*h)-exp(L2*h)-exp(-L2*h)*exp(L1*h)+exp(-L2*h)+exp(L1*h));
%a3 = exp(L2*h)*(1-exp(-c*dt))/((1-exp(-L1*h))*(1-exp(L2*h)));
a3 = (1-exp(-c*dt))*(exp(-L2*h)-exp(-L1*h))/(exp(L1*h)*exp(-L2*h)-exp(L1*h)-exp(-L2*h)-exp(L2*h)*exp(-L1*h)+exp(L2*h)+exp(-L1*h));
a2 = exp(-c*dt)-a1-a3;

T = fix((tf-ti)/dt);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);

%initial condition
u0 = zeros(N,1);
%exact solution
uex = zeros(N,1);
u = zeros(N,1);

for t=1:T
    for i=2:N-1
        u(i) = a1*u0(i-1)+a2*u0(i)+a3*u0(i+1)+(1-exp(-c*dt))*f(x(i),t*dt)/c;
    end
    %boundary conditions
    u(1) = exp(2*(x(1)-1)/Eps)*sin(2*t*dt) + exp(x(1)-1)*sin(t*dt);
    u(N) = exp(2*(x(N)-1)/Eps)*sin(2*t*dt) + exp(x(N)-1)*sin(t*dt);
    u0 = u;
end


%exact solution
for i=1:N
    uex(i,1) = exp(2*(x(i)-1)/Eps)*sin(2*tf) + exp(x(i)-1)*sin(tf);
end

if p==2
    plot(x,u,'*')
    hold on
    plot(x,uex)
    xlabel('x')
    ylabel('u')
    legend('approximate solution','exact solution')
end

err(p) = max(abs(u-uex));
h = h/2;

end

for i=1:p-1
    order(i) = log(err(i)/err(i+1))/log(2);
end
   

err
order