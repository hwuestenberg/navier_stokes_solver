%######################################################
% Apply a Successive over-relaxation method or Preconditioned Conjugate 
% Gradient method to compute pressure values for the next timestep
%
%
% LastLast: Added a function to correctly caluclate pressure terms depending on
% a cell's position (evaluate_p_terms)
%
%
% Last: Changed reshaping into vector AND reshaping back into a matrix.
% Columns and rows were exchanged in the function. Hence, one could not
% converge on rectangular grids (14/07/2018)
%
%
% Created on 25/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [P,it,res,t_solver] = Direct_Solver(P,RHS,K2D,imax,jmax)


        
%--------------------------- Direct Solver ------------------------------%
        
    % Reshape P matrix into a column vector of stacked columns from rows of P, similar for RHS
    rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);


    % Solve the linear system using a direct solver
    t_solver = tic;
    
    p = K2D\rhs;
    
    t_solver = toc(t_solver);
    res = norm(K2D*p-rhs,2);
    resvec = res;
    it=1;   % Unused because direct solver


    % Write pressure for next time step back into the matrix P
    P(2:imax,2:jmax) = reshape(p(:),[(imax-1) (jmax-1)]);

    % Free memory
    clear p rhs;
    
    
end
	
	