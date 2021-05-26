%######################################################
% Set boundary conditions at the edge of the computational grid
% Last: Added individual parameters for all boundaries (28/06/18)
%
% Created on 24/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [U,V,P,F,G] = boundary_values(U,V,P,imax,jmax,wl,wr,wt,wb)
%----------------------------- Velocity ----------------------------------
    % Left boundary (1,2:jmax)
    switch(wl)
        case 1 % Free-slip
            U(1,2:jmax) = 0;
            V(1,2:jmax) = V(2,2:jmax);
        case 2 % No-slip
            U(1,1:jmax) = 0;
            V(1,2:jmax) = -V(2,2:jmax);
        case 3 % Inflow
            U(1,2:jmax) = U(2,2:jmax);
            V(1,2:jmax) = V(2,2:jmax);
    end
    
    % Right boundary (imax+1,2:jmax)
    switch(wr)
        case 1 % Free-slip
            U(imax,2:jmax) = 0;
            V(imax+1,2:jmax) = V(imax,2:jmax);
        case 2 % No-slip
            U(imax,1:jmax) = 0;
            V(imax+1,2:jmax) = -V(imax,2:jmax);
        case 3 % Outflow
            U(imax,2:jmax) = U(imax-1,2:jmax);
            V(imax+1,2:jmax) = V(imax,2:jmax);
    end
    
    % Top boundary (2:imax,jmax+1)
    switch(wt)
        case 1 % Free-slip
            U(2:imax,jmax+1) = U(2:imax,jmax);
            V(2:imax,jmax) = 0;
        case 2 % No-slip
            U(2:imax,jmax+1) = -U(2:imax,jmax);
            V(jmax,1:imax) = 0;
        case 3 % Outflow
            U(2:imax,jmax+1) = U(2:imax,jmax);
            V(2:imax,jmax) = V(2:imax,jmax-1);
    end
    
    % Bottom boundary (2:imax,1)
    switch(wb)
        case 1 % Free-slip
            U(2:imax,1) = U(2:imax,2);
            V(2:imax,1) = 0;
        case 2 % No-slip
            U(2:imax,1) = -U(2:imax,2);
            V(1:imax,1) = 0;
        case 3 % Outflow
            U(2:imax,1) = U(2:imax,2);
            V(2:imax,1) = V(2:imax,2);
    end
    
%--------------------------- F and G terms -------------------------------    
    F = zeros(imax,jmax);
    G = zeros(imax,jmax);
    
    % Boundary conditions for Free-/No-slip and In-/Outflow
    F(1,2:jmax) = U(1,2:jmax);          % Left
    F(imax,2:jmax) = U(imax,2:jmax);    % Right
    G(2:imax,1) = V(2:imax,1);          % Bottom
    G(2:imax,jmax) = V(2:imax,jmax);    % Top
    
%----------------------------- Pressure ----------------------------------
    % Boundary conditions for Free-/No-slip and In-/Outflow
    P(1,2:jmax) = P(2,2:jmax);          % Left
    P(imax+1,2:jmax) = P(imax,2:jmax);  % Right
    P(2:imax,jmax+1) = P(2:imax,jmax);  % Top
    P(2:imax,1) = P(2:imax,2);          % Bottom
end



