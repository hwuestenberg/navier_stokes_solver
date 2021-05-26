%  SOR method to approximate x of Ax=y.
%
%  Author: Martin Krosche
%          Institute of Scientific Computing
%          Technische Universitaet Braunschweig
%          Braunschweig, Germany
%
%  Copyright (c) 2008 Martin Krosche. All rights reserved. No warranty. No
%  liability.
%
%  Original method slightly changed to measure iterations and the residium 
%  during every iteration for comparison.
%  by Henrik Wüstenberg with Matlab R2017b
%
%  [X, ERR] = SOR( A, Y, X0, OMEGA, ITER ) returns X by performing ITER
%  iterations and taking starting vector X0. ERR is the L2-norm of the residuum.
function [x, i, err, resvec] = sor_algorithm( A, y, x0, omega, itermax, tol )
    
    % Intitialise vector of residuum norms
    resvec = zeros(itermax,1);
    
    %  get diagonal in sparse format.
    n    = length( A );
    D    = spdiags( diag( A ), 0, n, n );

    %  get lower and upper triangular.
    L = tril( A, -1 );
    U = triu( A, 1 );

    B = inv( D + omega*L );
    M = B * ( (1-omega)*D - omega*U );
    g = B * omega*y;

    %  iterate...
    x = x0;
    for i=1:itermax
        err = norm( ( A*x-y ), 2 );
        resvec(i) = err;
        if err < tol
            break;
        end
        x = M*x + g;
    end
