%######################################################
% Initialize the computational grid for solving a given fluid problem. The
% initial parameters are equally set for each element on the grid.
%
% Created on 24/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [U,V,P] = initgrid(imax,jmax,U_I,V_I,P_I)
    % Set initial values on all grid positions
    U = ones(imax+1,jmax+1)*U_I;
    V = ones(imax+1,jmax+1)*V_I;
    P = ones(imax+1,jmax+1)*P_I;
end

