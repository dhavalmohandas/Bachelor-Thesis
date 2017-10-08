clear all
close all
clc
format long
xi = 1;
xf = 3;
h = 1/(2^6);
ti = 0;
tf = 0.5;
dt = 1/(2^10);
T = fix((tf-ti)/dt);

for c=1:4
N = fix((xf-xi)/h)+1;
x = linspace(xi,xf,N);

u0 = zeros(N,1);
u = zeros(N,1);

%initial condition
for j=1:N
    u0(j) = exp(-50*x(j)^2)*exp( 1i*sin(x(j)) );
end

a = pi + x;
k = cos(x);
A1 = zeros(N-2);
A2 = zeros(N-2);


for m=1:T
    
    for j=2:N-1
        a1(j-1) = sin(k(j)*(a(j)*dt-h)/2)/sin(k(j)*(a(j)*dt+h)/2);
        a2(j-1) = 1;
        A1(j-1,j-1) = 1;
        A2(j-1,j-1) = -a1(j-1);
        if j>2
            A1(j-1,j-2) = -a1(j-1);
            A2(j-1,j-2) = a2(j-1);
        end
    end
    u(1) = exp(-50*(a(1)*exp(-m*dt)-pi)^2)*exp(1i*sin(a(1)*exp(-m*dt)-pi));
    u(N) = exp(-50*(a(N)*exp(-m*dt)-pi)^2)*exp(1i*sin(a(N)*exp(-m*dt)-pi));
    %u(1) = exp(-50*( x(1)-a(1)*tf)^2)*exp( 1i*sin( x(1)-a(1)*tf ) );
    %u(N) = exp(-50*( x(N)-a(N)*tf)^2)*exp( 1i*sin( x(N)-a(N)*tf ) );
    R = A2*u0(2:N-1);
    R(1) = R(1) + a1(1)*u(1)+a2(1)*u0(1);
    u(2:N-1) = A1\R;
    u0 = u;
end


uex = zeros(N,1); %exact solution

for j=1:N
    uex(j,1) = exp(-50*( a(j)*exp(-tf)-pi )^2)*exp( 1i*sin( a(j)*exp(-tf)-pi ) );
    %uex(j,1) = exp(-50*( x(j)-a(j)*tf)^2)*exp( 1i*sin( x(j)-a(j)*tf ) );
end

if c==2
    figure(1)
    plot(x,real(u),'*')
    hold on
    plot(x,real(uex))
    ylabel('real(u)')
    xlabel('x')
    title('real part')
    legend('approx solution', 'exact solution')
    
    figure(2)
    plot(x,imag(u),'*')
    hold on
    plot(x,imag(uex))
    ylabel('imag(u)')
    xlabel('x')
    title('imaginary part')
    legend('approx solution', 'exact solution')
end



err(c) = max(abs(u-uex));
h = h/2;
end

for j=1:c-1
    order(j) = log(err(j)/err(j+1))/log(2);
end





order
err
