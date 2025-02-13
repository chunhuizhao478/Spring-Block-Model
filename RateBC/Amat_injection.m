%Generate stencil
function [A_mat] = Amat_injection(x,D_c,Lx) %D_c = L / v_o / (phi * ct * mu )
   n = length(x);
   dx = x(2)-x(1);
   %get left,right k value for each node
   k_list = zeros(length(x),2);
   for i = 1:n
       k_left = k(x(i)-0.5);
       k_right = k(x(i)+0.5);
       k_list(i,1) = k_left;
       k_list(i,2) = k_right;
   end
   %make value list for low,mid,upp
   low_tridiag = zeros(n,1);
   mid_tridiag = zeros(n,1);
   upp_tridiag = zeros(n,1);
   for j = 1:n
       low_tridiag(j) = D_c * k_list(j,1) / (dx^2 * Lx^2);                
       mid_tridiag(j) = - D_c * (k_list(j,1) + k_list(j,2)) / (dx^2 * Lx^2);
       upp_tridiag(j) = D_c * k_list(j,2) / (dx^2 * Lx^2);                     
   end
   %construct A mat
   %assume Dirichlet condition at both ends, so the size of mat is (n-2)by(n-2)
   A_mat = (diag(low_tridiag(3:n-1),-1) + diag(mid_tridiag(2:n-1)) + diag(upp_tridiag(2:n-2),1));
   %Change to flow rate boundary condition
   %%calculate kP_n+1'' kP_n-1''
   %Multiply Nondimensionalize factor
%    A_mat(ceil(end/2),:) = 0;
%    A_mat(ceil(end/2),ceil(end/2)-1) =  - k(0) / (mu * dx) * factor;
%    A_mat(ceil(end/2),ceil(end/2)+0) =  0;
%    A_mat(ceil(end/2),ceil(end/2)+1) =  k(0) / (mu * dx) * factor;
end