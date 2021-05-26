%######################################################
% Set boundary conditions at the edge of the computational grid
% Last: Created grid-wise calculation
%
% Created on 12/06/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [Psi,Zeta] = psi_zeta(U,V,imax,jmax,delx,dely)
    Psi = zeros(imax,jmax);
    Zeta = zeros(imax,jmax);
    
    % Predefinition of Psi at the lower boundary (at i=0 instead of i=1)
    Psi(2:imax,1) = 0;
    
    % Compute Psi over the grid
    for i=2:imax
        for j=2:jmax
            Psi(i,j) = U(i,j)*dely+Psi(i,j-1);
            Zeta(i,j) = (U(i,j)-U(i,j-1))/dely - (V(i,j)-V(i,j-1))/delx;
        end
    end
end



