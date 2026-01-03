function [x,Dx,Dxx] = diffcheb(n,xspan)

x = -cos( (0:n)'*pi/n);
Dx = zeros(n+1);
c = [2; ones(n-1,1); 2];
i = (0:n)';

for j = 0:n
    num = c(i+1).*(-1).^(i+j);
    den = c(j+1)*(x-x(j+1));
    Dx(:,j+1) = num./den;
end

Dx(isinf(Dx)) = 0;
Dx = Dx - diag(sum(Dx,2));

a = xspan(1); b = xspan(2);
x = a + (b-a)*(x+1)/2;
Dx = 2*Dx/(b-a);

Dxx = Dx^2;

end