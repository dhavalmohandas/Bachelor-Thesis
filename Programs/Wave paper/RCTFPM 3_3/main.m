clear all
close all

format long
xi = -0.5;
xf = 0.5;
ti = 0;
tf = 0.4;
dt = 1/2^6;

h = 1/2^7;
epsilon = 0.02;

for z=1:4


    N = fix((xf-xi)/h)+1;
    T = fix((tf-ti)/dt);
    x = linspace(xi,xf,N);

    u0 = zeros(N,1);
    u = zeros(N,1);

    %initial condition
    for j=1:N
        u0(j) = exp(-10*(2*x(j)+cos(x(j)))^2)*exp(1i*sin(2*x(j)+cos(x(j)) )/epsilon);
    end



    for j=1:N
        a(j) = (2*x(j)+cos(x(j)))/(2-sin(x(j)));
    end

    A1 = zeros(N-2);
    A2 = zeros(N-2);
    
    
    for m=1:T

        for j=2:N-1
            b = (a(j-1)*x(j)-a(j)*x(j-1))/h;
            c = (a(j)-a(j-1))/h;
            y = ((b + c*x(j))*exp(-m*dt)-b)/c;
            k = cos(y)*exp(-c*m*dt)/epsilon;
            a1(j-1) = sin(k*(a(j)*dt-h)/2)/sin(k*(a(j)*dt+h)/2);
            a2(j-1) = 1;
            A1(j-1,j-1) = 1;
            A2(j-1,j-1) = -a1(j-1);
            if j>2
                A1(j-1,j-2) = -a1(j-1);
                A2(j-1,j-2) = a2(j-1);
            end        
        end
        u(1) = exp(-10*((2*x(1)+cos(x(1)))*exp(-m*dt))^2)*exp(1i*sin((2*x(1)+cos(x(1)))*exp(-m*dt))/epsilon);
        u(N) = exp(-10*((2*x(N)+cos(x(N)))*exp(-m*dt))^2)*exp(1i*sin((2*x(N)+cos(x(N)))*exp(-m*dt))/epsilon);
        R = A2*u0(2:N-1);
        R(1) = R(1) + a1(1)*u(1)+a2(1)*u0(1);
        u(2:N-1) = A1\R;
        u0 = u;

    end

    usol = zeros(N,1);
    %exact solution
    for j=1:N
        usol(j) = exp(-10*((2*x(j)+cos(x(j)))*exp(-tf))^2)*exp(1i*sin((2*x(j)+cos(x(j)))*exp(-tf))/epsilon);
    end

    err(z) = max(max(abs(u-usol)));
    h = h/2;
    dt = dt/2;
end

for j=1:z-1
    order(j) = log(err(j)/err(j+1))/log(2);
end


figure(1)
plot(x,real(usol))
hold on
plot(x,real(u))
title('real part')
legend('exact solution','approx solution')

figure(2)
plot(x,imag(usol))
hold on
plot(x,imag(u))
title('imaginary part')
legend('exact solution','approx solution')

err
order
