clear all
close all

format long
xi = 0;
xf = 1;
ti = 0;
tf = 0.25;
dt = 1/2^6;

h = 1/2^7;
epsilon = 0.01;

for c=1:4


    N = fix((xf-xi)/h)+1;
    T = fix((tf-ti)/dt);
    x = linspace(xi,xf,N);

    u0 = zeros(N,1);
    u = zeros(N,1);

    %initial condition
    for j=1:N
        u0(j) = exp(-200*x(j)^2)*exp(1i*sin(x(j))/epsilon);
    end



    a = 0.5 + x;

    A1 = zeros(N-2);
    A2 = zeros(N-2);
    
    for m=1:T

        for j=2:N-1
            y(j) = (0.5 + x(j))*exp(-m*dt)-0.5;
            k(j) = cos(y(j))*exp(-m*dt)/epsilon;
            a1(j-1) = sin(k(j)*(a(j)*dt-h)/2)/sin(k(j)*(a(j)*dt+h)/2);
            a2(j-1) = 1;
            A1(j-1,j-1) = 1;
            A2(j-1,j-1) = -a1(j-1);
            if j>2
                A1(j-1,j-2) = -a1(j-1);
                A2(j-1,j-2) = a2(j-1);
            end        
        end
        %u(1) = exp(-200*(x(1)-a(1)*m*dt)^2)*exp(1i*sin(x(1)-a(1)*m*dt)/epsilon); %boundary condition
        %u(N) = exp(-200*(x(N)-a(N)*m*dt)^2)*exp(1i*sin(x(N)-a(N)*m*dt)/epsilon);
        u(1) = exp(-200*(a(1)*exp(-m*dt)-0.5)^2)*exp(1i*sin(a(1)*exp(-m*dt)-0.5)/epsilon);
        u(N) = exp(-200*(a(N)*exp(-m*dt)-0.5)^2)*exp(1i*sin(a(N)*exp(-m*dt)-0.5)/epsilon);
        R = A2*u0(2:N-1);
        R(1) = R(1) + a1(1)*u(1)+a2(1)*u0(1);
        u(2:N-1) = A1\R;
        u0 = u;

    end

    usol = zeros(N,1);
    %exact solution
    for j=1:N
        %usol(j) = exp(-200*(x(j)-a(j)*tf)^2)*exp(1i*sin(x(j)-a(j)*tf)/epsilon);
        usol(j) = exp(-200*(a(j)*exp(-tf)-0.5)^2)*exp(1i*sin(a(j)*exp(-tf)-0.5)/epsilon);
    end

    err(c) = max(max(abs(u-usol)));
    h = h/2;
    dt = dt/2;
end

for j=1:c-1
    order(j) = log(err(j)/err(j+1))/log(2);
end


figure(1)
plot(x,real(usol))
hold on
plot(x,real(u),'*')
xlabel('x')
ylabel('real(u)')
title('real part')
legend('exact solution','approx solution')

figure(2)
plot(x,imag(usol))
hold on
plot(x,imag(u),'*')
xlabel('x')
ylabel('imag(u)')
title('imaginary part')
legend('exact solution','approx solution')
err
order
