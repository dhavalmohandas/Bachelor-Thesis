clear all
close all
format long
xi = 0;
xf = 1;
ti = 0;
tf = 1;
h = 1/2^5;
dt = 2^-6;
T = fix((tf-ti)/dt);
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);
eps = 0.01;

%intial condition
u0 = zeros(N,1);

u = zeros(N,1);
uex = zeros(N,1);
A1 = zeros(N-2);
A2 = zeros(N-2);

c = 0;

for i=1:N
    if x(i)<=0.5
        b(i) = 1;
    else
        b(i) = -1;
    end
end

for t=1:T
    %boundary conditions
    u(1) = (t*dt)^2;
    u(N) = 0;
    for j=2:N-1
        a1 = (exp(b(j)*h/eps-c*dt/2)*b(j)*dt)/(2*h*(exp(b(j)*h/eps)-1));
        a3 = (exp(-c*dt/2)*b(j)*dt)/(2*h*(exp(b(j)*h/eps)-1));
        a2 = exp(-c*dt/2)-a1-a3;
        b1 = -(exp(b(j)*h/eps+c*dt/2)*b(j)*dt)/(2*h*(exp(b(j)*h/eps)-1));
        b3 = -(exp(c*dt/2)*b(j)*dt)/(2*h*(exp(b(j)*h/eps)-1));
        b2 = exp(c*dt/2)-b1-b3;
        A1(j-1,j-1) = a2;
        A2(j-1,j-1) = b2;
        if j>2
            A1(j-1,j-2) = a1;
            A2(j-1,j-2) = b1;
        end
        if j<N-2
            A1(j-1,j) = a3;
            A2(j-1,j) = b3;
        end
    end
    R = A1*u0(2:N-1);
    R(1) = R(1) + a1*u0(1)-b1*u(1);
    R(N-2) = R(N-2) + a3*u(N)-b3*u(N);
    u(2:N-1) = A2\R;
    
    u0 = u;
end


plot(x,u)