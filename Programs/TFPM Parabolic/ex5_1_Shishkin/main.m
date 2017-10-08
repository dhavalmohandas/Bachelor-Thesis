clear all
close all
format long
clc
global Eps 
xi = 0;
xf = 1;
ti = 0;
tf = 1;
b = 1;
Eps = 0.01;  
N = 32;
dt = 1/4096;
T = fix((tf-ti)/dt);

%for p=1:3
    
h = (xf-xi)/N;
%shishkin mesh
sigma = max(0.5,1-2*Eps*log(N)/b);
x1 = linspace(xi,sigma,fix(N/2)+1);
x2 = linspace(sigma,xf,fix(N/2)+1);
x = cat(2,x1,x2(2:fix(N/2)+1));
A1 = zeros(N-1);

%initial condition
v0 = zeros(N+1,1);
v = zeros(N+1,1);
u = zeros(N+1,1);

for t=1:T
    v(1) = exp((x(1)-1)/Eps)*sin(2*t*dt)+exp((x(1)-1)/2)*sin(t*dt)-f(x(1),t*dt)*t*dt;
    v(N+1) = exp((x(N+1)-1)/Eps)*sin(2*t*dt)+exp((x(N+1)-1)/2)*sin(t*dt)-f(x(N+1),t*dt)*t*dt;
    for j=2:N
        h = x(j)-x(j-1);
        h1 = x(j+1)-x(j);
        a1 = b*dt*(exp(b*h1/Eps)-1)/(h1*(1-exp(-b*h/Eps))-h*(exp(b*h1/Eps)-1));
        a3 = b*dt*(1-exp(-b*h/Eps))/(h1*(1-exp(-b*h/Eps))-h*(exp(b*h1/Eps)-1));
        a2 = 1-a1-a3;
        A(j-1,j-1) = a2;
        if j>2
            A(j-1,j-2) = a1;
        end
        if j<N
            A(j-1,j) = a3;
        end
    end
    a1 = b*dt*(exp(b*(x(3)-x(2))/Eps)-1)/((x(3)-x(2))*(1-exp(-b*(x(2)-x(1))/Eps))-(x(2)-x(1))*(exp(b*(x(3)-x(2))/Eps)-1));
    a3 = b*dt*(1-exp(-b*(x(N)-x(N-1))/Eps))/((x(N+1)-x(N))*(1-exp(-b*(x(N)-x(N-1))/Eps))-(x(N)-x(N-1))*(exp(b*(x(N+1)-x(N))/Eps)-1));
    
    R = v0(2:N);
    R(1) = R(1)-a1*v(1);
    R(N-1) = R(N-1)-a3*v(N+1);
    v(2:N) = A\R;
    v0 = v;
end

for i=1:N+1
    u(i,1) = v(i)+f(x(i),tf)*tf;
end

%exact solution;
for i=1:N+1
    uex(i,1) = exp((x(i)-1)/Eps)*sin(2*tf)+exp((x(i)-1)/2)*sin(tf);
end



plot(x,uex);
hold on
plot(x,u)
xlabel('x')
ylabel('y')
legend('exact solution','approximate solution')

