clear all
close all
clc
format long

xi = 0;
xf = 2;
ti = 0;
tf = 1;
h = 1/2^7;
dt = 1/2^8;
N = fix((xf-xi)/h)+1;
T = fix((tf-ti)/dt);
x = linspace(xi,xf,N);
lambda = dt/h;
u0 = zeros(1,N);
u = zeros(1,N);

%initial condition
for i=1:N
    if x(i)>=0.4 && x(i)<=0.6
        u0(i) = 1;
    else
        u0(i) = 0;
    end
end

for i=1:1
    a(1) = (sin((u0(1)+u0(2))*dt/4)*sin(((u0(2) + u0(1))*dt/2-h)/2))/(sin(h/2)*sin(h));
    b(1) = (sin((u0(1)+u0(2))*dt/4)*sin(((u0(2) + u0(1))*dt/2+h)/2))/(sin(h/2)*sin(h));
    for j=2:N-1
        a(j) = (sin((u0(j)+u0(j+1))*dt/4)*sin(((u0(j+1) + u0(j))*dt/2-h)/2))/(sin(h/2)*sin(h));
        b(j) = (sin((u0(j)+u0(j+1))*dt/4)*sin(((u0(j+1) + u0(j))*dt/2+h)/2))/(sin(h/2)*sin(h));
        theta = (u0(j)-u0(j-1))/(u0(j+1)-u0(j));
        phi = (theta + abs(theta))/(1+abs(theta));
        g1 = u0(j)^2/2 + phi*(-(a(j)*u0(j+1)-b(j)*u0(j))/lambda -(u0(j+1)+u0(j))*u0(j)/2 );
        g2 = u0(j-1)^2/2 + phi*(-(a(j-1)*u0(j)-b(j-1)*u0(j-1))/lambda -(u0(j-1)+u0(j))*u0(j-1)/2 );
        u(j) = u0(j)-lambda*(g1-g2);
    end
    u0 = u;
end

usol = zeros(1,N);

for i=1:N
    if x(i)<0.4;
        usol(i) = 0;
    elseif x(i)< 0.4 + 2*sqrt(tf)/sqrt(10)
        usol(i) = (x(i)-0.4)/tf;
    else
        usol(i) = 0;
    end
end

figure(1)
plot(x,usol)
figure(2)
plot(x,u)
