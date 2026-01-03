function [x,Dx,Dxx] = diffmat2(n,xspan)

x  = linspace(xspan(1),xspan(end),n);
dx = x(2) - x(1);

% First Derivative Matrix;
Dx = diag(1/2*ones(1,n-1),1) - ...
     diag(1/2*ones(1,n-1),-1); 
Dx(1,:) = [-3/2 2 -1/2 zeros(1,n-3)];
Dx(end,:) = [zeros(1,n-3) 3/2 -2 1/2];
Dx = Dx/dx;

% Second Derivative Matrix;
Dxx = diag(ones(1,n-1),1) + ...
      diag(ones(1,n-1),-1) - ...
      diag(2*ones(1,n));
Dxx(1,:) = [2 -5 4 -1 zeros(1,n-4)];
Dxx(end,:) = [zeros(1,n-4) -1 4 -5 2];
Dxx = Dxx/dx/dx;

Dx  = sparse(Dx);
Dxx = sparse(Dxx);


end
