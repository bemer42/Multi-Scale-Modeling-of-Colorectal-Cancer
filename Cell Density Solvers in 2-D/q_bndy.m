function Boundary = q_bndy(N_x,N_y)

% Define Entire Boundary:
% Boundary.isbndy = true(N_x,N_y);
% Boundary.isbndy(2:N_x-1,2:N_y-1) = false;

% Define the top, bottom, left, and right boundary:
Boundary.istop = true(N_x,N_y);
Boundary.istop(2:N_x,1:N_y) = false;
Boundary.istop(1,1) = false;
Boundary.istop(1,N_y) = false;

Boundary.isbottom = true(N_x,N_y);
Boundary.isbottom(1:N_x-1,1:N_y) = false;
Boundary.isbottom(N_x,1) = false;
Boundary.isbottom(N_x,N_y) = false;

Boundary.isleft = true(N_x,N_y);
Boundary.isleft(1:N_x,2:N_y) = false;
Boundary.isleft(1,1) = false;
Boundary.isleft(N_x,1) = false;

Boundary.isright = true(N_x,N_y);
Boundary.isright(1:N_x,1:N_y-1) = false;
Boundary.isright(1,N_y) = false;
Boundary.isright(N_x,N_y) = false;

% Define the TL, TR, BL, and BR corner:
Boundary.isTL = false(N_x,N_y);
Boundary.isTL(1,1) = true;

Boundary.isTR = false(N_x,N_y);
Boundary.isTR(1,N_y) = true;

Boundary.isBL = false(N_x,N_y);
Boundary.isBL(N_x,1) = true;

Boundary.isBR = false(N_x,N_y);
Boundary.isBR(N_x,N_y) = true;

% Full Boundary:
Boundary.isbndy = Boundary.istop + Boundary.isbottom + Boundary.isleft + ...
                  Boundary.isright + Boundary.isTL + Boundary.isTR + ...
                  Boundary.isBL + Boundary.isBR;

Boundary.vec = @(U) U(:);
Boundary.unvec = @(u) reshape(u,N_x,N_y);

end