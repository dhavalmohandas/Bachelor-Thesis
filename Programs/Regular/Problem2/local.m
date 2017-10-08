function coeff = local(p,q) 
xi = 0;
xf = 2;
yi = 0;
yf = 1;
M = 5;   % Number of coefficients/terms in scheme
N = 5;   % Number of points in the uniform mesh
x1 = linspace(xi,xf,2*N);
y1 = linspace(yi,yf,N);
h = (xf-xi)/(N-1);

syms a b


basis = [1,a,b,a^2-b^2];

for i=1:M
    
    if i<5
        u = basis(i);
        A(i,1) = subs(u,a,x1(p-1));
        A(i,1) = subs(A(i,1),b,y1(q));

        A(i,2) = subs(u,a,x1(p+1));
        A(i,2) = subs(A(i,2),b,y1(q));

        A(i,3) = subs(u,a,x1(p));
        A(i,3) = subs(A(i,3),b,y1(q));

        A(i,4) = subs(u,a,x1(p));
        A(i,4) = subs(A(i,4),b,y1(q-1));
    
    
        A(i,5) = subs(u,a,x1(p));
        A(i,5) = subs(A(i,5),b,y1(q+1));
    else
        A(i,1) = 0;
        A(i,2) = 0;
        %A(i,3) = -4/h^2;
        A(i,3) = 1;
        A(i,4) = 0;
        A(i,5) = 0;
    end
    
end

%v(1) = f(x1(p-1,y1(q)));

v = f(x1(p),y1(q))*zeros(5,1);
v(5) = -4/h^2;

%v = zeros(5,1);
%v(5) = -4/h^2;

coeff = double(A\v);
end