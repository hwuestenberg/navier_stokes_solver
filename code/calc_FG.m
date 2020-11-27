%######################################################
% Calculate F and G. These terms are used for solving the Poisson equation
% and computing velocities for the next time step in this simulation.
% Last: Matrices F and G are created in function 'boundary_values'
%
% Created on 24/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [F,G] = calc_FG(U,V,F,G,imax,jmax,delt,delx,dely,GX,GY,alpha,nu)
    % Compute F and G (Compute column by column (j-wise))
    for i = 2:imax-1
        for j = 2:jmax
            d2u = 1/delx^2*( U(i+1,j)-2*U(i,j)+U(i-1,j) ) + 1/dely^2*( U(i,j+1)-2*U(i,j)+U(i,j-1) );
            duu_dx = 1/delx*( ( ( U(i,j)+U(i+1,j) )/2 )^2 - ( ( U(i-1,j)+U(i,j) )/2 )^2 ) ...
                + alpha/delx*( ( abs(U(i,j)+U(i+1,j)) )/2 * ( U(i,j)-U(i+1,j) )/2 - ( abs(U(i-1,j)+U(i,j)) )/2 * ( U(i-1,j)-U(i,j) )/2 );
            duv_dy = 1/dely*( ( V(i,j)+V(i+1,j) )/2 * ( U(i,j)+U(i,j+1) )/2 - ( V(i,j-1)+V(i+1,j-1) )/2 * ( U(i,j-1)+U(i,j) )/2 ) ...
                + alpha/dely*( ( abs(V(i,j)+V(i+1,j)) )/2 * ( U(i,j)-U(i,j+1) )/2 - ( abs(V(i,j-1)+V(i+1,j-1)) )/2 * ( U(i,j-1)-U(i,j) )/2 );
            Fbracket = nu*(d2u) - duu_dx - duv_dy + GX;
            F(i,j) = U(i,j) + delt*Fbracket;
		end
	end
	
	for i = 2:imax
		for j = 2:jmax-1
            d2v = 1/delx^2*( V(i+1,j)-2*V(i,j)+V(i-1,j) ) + 1/dely^2*( V(i,j+1)-2*V(i,j)+V(i,j-1) );
            dvv_dy = 1/dely*( ( ( V(i,j)+V(i,j+1) )/2 )^2 - ( ( V(i,j-1)+V(i,j) )/2 )^2 ) ...
                + alpha/dely*( ( abs(V(i,j)+V(i,j+1)) )/2 * ( V(i,j)-V(i,j+1) )/2 - ( abs(V(i,j-1)+V(i,j)) )/2 * ( V(i,j-1)-V(i,j) )/2 );
            duv_dx = 1/delx*( ( U(i,j)+U(i,j+1) )/2 * ( V(i,j)+V(i+1,j) )/2 - ( U(i-1,j)+U(i-1,j+1) )/2 * ( V(i-1,j)+V(i,j) )/2 ) ...
                + alpha/delx*( ( abs(U(i,j)+U(i,j+1)) )/2 * ( V(i,j)-V(i+1,j) )/2 - ( abs(U(i-1,j)+U(i-1,j+1)) )/2 * ( V(i-1,j)-V(i,j) )/2 );
            Gbracket = nu*(d2v) - dvv_dy - duv_dx + GY;
            G(i,j) = V(i,j) + delt*Gbracket;
        end
    end
    
end



