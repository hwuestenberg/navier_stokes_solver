%######################################################
% Create a Block tridiagonal matrix for solving a linear system of a
% two-dimensional problem in function 'SOR'
%
% Created on 12/06/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function K2D = init_linSys(imax,jmax,dely,delx)
    % Set up a linear system of equations
    vy = sparse([2 -1 zeros(1,jmax-3)]);    % -3 for both numbers and the +1 incorporated in jmax
    K1Dy = (1/dely^2)*toeplitz(vy);         % Create Toeplitz matrix for one dimension
    Iy = speye(jmax-1);                     % Identity matrix for block tridiagonal matrix
%     fprintf('Maximum eigenvalue of K1Dy matrix is %f\n',max(eig(K1Dy)));
    
    vx = sparse([2 -1 zeros(1,imax-3)]);
    K1Dx = (1/delx^2)*toeplitz(vx);
    Ix = speye(imax-1);
%     fprintf('Maximum eigenvalue of K1Dx matrix is %f\n',max(eig(K1Dx)));

    K2D = kron(K1Dy,Ix)+kron(Iy,K1Dx);  % Create block tridiagonal matrix for two dimensions
    clear K1Dx K1Dy Ix Iy vx vy         % Free memory
end
