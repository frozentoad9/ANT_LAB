h1 = 0.25;
x1=(0:h1:1)';
h2 = 0.1;
x2=(0:h2:1)';
y1 = ff(h1);
y2 = ff(h2);
plot(y1,x1,'b');
hold on
plot(y2,x2,'r');

legend('plot for h=0.25','plot for h=0.1',Location='best')
grid on

function y = ff(h)
    x=(0:h:1)';
    m=length(x);
    y=zeros(m,1);
    y(m,1)=2;
    
    
    A=ones(m,1);
    B=zeros(m,1);
    C=x.*(-1);
    D=ones(m,1).*(2);
    
    
    
    a=(A./(h^2))-(B./(2*h));
    b=((A./(h^2)).*(-2))+C;
    c=(A./(h^2))+(B./(2*h));
    d=D;
    
    
    
    alpha_0=1.0;
    alpha_1=1.0;
    gamma_1=2.0;
    
    
    
    Input=zeros(m-1,m-1); 
    
    
    
    Input(1,1)=b(1,1)+((alpha_0/alpha_1))*2*h * a(1,1);
    Input(1,2)=c(1,1)+1*a(1,1);
    
    
    
    d(1,1)=D(1,1)+((a(1,1)*2.0*h*gamma_1)/alpha_1);
    
    
    
    i=2;
    while(i<m-1)
    Input(i,i)=b(i,1);
    Input(i,i-1)=a(i,1);
    Input(i,i+1)=c(i,1);
    i=i+1;
    end
    
    
    
    Input(m-1,m-1) = b(m-1,1);
    Input(m-1,m-2)=a(m-1,1);
    d(m-1,1)=D(m-1,1)-(c(m-1,1)*y(m,1));
    
    
    
    i=1;
    RHS=zeros(m-1,1);
    while(i<=m-1)
    RHS(i,1)=d(i,1);
    i=i+1;
    end
   
    y_temp=Gauss_elem(Input,RHS,m-1);
    
    
    
    i=1;
    while(i<m)
    y(i,1)=y_temp(i,1);
    i=i+1;
    end
    
   
    fprintf("Output vector for h = %f is (y0 to yn): \n",h);
    display(y(:,1));
end

 function x = Gauss_elem(Input,d,n)
    A = [Input d];                                                  %Augmented Matrix
    x = zeros(n,1);                                             %variable matrix [x1 x2 ... xn] coulmn
    for i=1:n-1
        for j=i+1:n
            m = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - m*A(i,:);
        end
    end
    x(n) = A(n,n+1)/A(n,n);
    for i=n-1:-1:1
        summ = 0;
        for j=i+1:n
            summ = summ + A(i,j)*x(j,:);
        end
        x(i,:) = (A(i,n+1) - summ)/A(i,i);
    end
 end

%  Solution is 
% Output vector for h = 0.250000 is (y0 to yn): 
%   -14.9232
%   -10.6299
%    -6.3777
%    -2.1998
%     2.0000
% 
% Output vector for h = 0.100000 is (y0 to yn): 
%   -14.1796
%   -12.5516
%   -10.9162
%    -9.2826
%    -7.6569
%    -6.0418
%    -4.4369
%    -2.8387
%    -1.2403
%     0.3682
%     2.0000




