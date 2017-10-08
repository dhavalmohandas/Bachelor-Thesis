function coeff = local(s,t)

xi = 0;
xf = 1;
yi = 0;
yf = 1;
epsilon = 1;

N = 6; % number of points in the mesh

M = 5; % number of terms in the scheme

x1 = linspace(xi,xf,N);
y1 = linspace(yi,yf,N);

h = (xf-xi)/(N-1);

p = 1;
q = 0;
b = 1;

global mu

mu = sqrt(b + (p^2+q^2)/(4*epsilon^2))/epsilon;

A = zeros(M);

for i=1:M
    if i<5
        
        A(i,1) = basis(i,x1(s+1),y1(t));        

        A(i,2) = basis(i,x1(s-1),y1(t));        

        A(i,3) = basis(i,x1(s),y1(t));        

        A(i,4) = basis(i,x1(s),y1(t+1));    
    
        A(i,5) = basis(i,x1(s),y1(t-1));
        
    else
        A(i,1) = 0;
        A(i,2) = 0;
        A(i,3) = 1;
        A(i,4) = 0;
        A(i,5) = 0;
    end
        
end


v = zeros(5,1);
v(5) = (cosh(0.5*mu*h)/sinh(0.5*mu*h))^2;


coeff = A\v;
end