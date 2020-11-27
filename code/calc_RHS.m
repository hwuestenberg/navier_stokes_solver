%######################################################
% Calculate the right hand side of the pressure equation based on the terms
% F and G
%
% Last: Removed a function to correctly calculate the rhs for each cell's
% position (causes "early" divergence)
%
% Created on 25/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function RHS = calc_RHS(F,G,imax,jmax,delt,delx,dely)
    RHS = zeros(imax,jmax);
    
    % Compute right hand side
    for i = 2:imax
        for j = 2:jmax
            f = F(i-1,j)-F(i,j);
            g = G(i,j-1)-G(i,j);
            RHS(i,j) = 1/delt*( 1/delx*( f ) + 1/dely*( g ) );
        end
    end
end



