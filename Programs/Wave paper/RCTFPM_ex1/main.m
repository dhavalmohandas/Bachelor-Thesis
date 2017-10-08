clear all
close all


format long
xi = 0;
xf = 2;
ti = 0;
tf = 1;
dt = 1/2^10;
T = (tf-ti)/(dt);
h = 1/2^6;

for c=1:4


    N = fix((xf-xi)/h)+1;

    x = linspace(xi,xf,N);

    u0 = zeros(N,1);
    u = zeros(N,1);

    %initial condition
    for j=1:N
        u0(j) = exp(-50*x(j)^2)*exp(1i*sin(x(j)));
    end



    k = cos(x);

    A1 = zeros(N-2);
    A2 = zeros(N-2);

    for m=1:T

        for j=2:N-1
            a1(j-1) = sin(k(j)*(dt-h)/2)/sin(k(j)*(dt+h)/2);
            a2(j-1) = 1;
            A1(j-1,j-1) = 1;
            A2(j-1,j-1) = -a1(j-1);
            if j>2
                A1(j-1,j-2) = -a1(j-1);
                A2(j-1,j-2) = a2(j-1);
            end        
        end
        u(1) = exp(-50*(x(1)-m*dt)^2)*exp(1i*sin(x(1)-m*dt)); %boundary condition
        u(N) = exp(-50*(x(N)-m*dt)^2)*exp(1i*sin(x(N)-m*dt));
        R = A2*u0(2:N-1);
        R(1) = R(1) + a1(1)*u(1)+a2(1)*u0(1);
        u(2:N-1) = A1\R;
        u0 = u;

    end

    usol = zeros(N,1);
    %exact solution
    for j=1:N
        usol(j) = exp(-50*(x(j)-1)^2)*exp(1i*sin(x(j)-1));
    end

    err(c) = max(max(abs(u-usol)));
    h = h/2;
    if c==2
        figure(1)
        plot(x,real(usol))
        hold on
        plot(x,real(u),'*')
        title('real part')
        xlabel('x')
        ylabel('real(u)')
        legend('exact solution','approx solution')

        figure(2)
        plot(x,imag(usol))
        hold on
        plot(x,imag(u),'*')
        title('imaginary part')
        xlabel('x')
        ylabel('imag(u)')
        legend('exact solution','approx solution')
    end 
    
end

for j=1:c-1
    order(j) = log(err(j)/err(j+1))/log(2);
end



err
order
