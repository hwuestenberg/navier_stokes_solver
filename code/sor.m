%######################################################
% Apply a Successive over-relaxation method or Preconditioned Conjugate 
% Gradient method to compute pressure values for the next timestep
% LastLast: Added a function to correctly caluclate pressure terms depending on
% a cell's position
% Last: Changed reshaping into vector AND reshaping back into a matrix.
% Columns and rows were exchanged in the function. Hence, one could not
% converge on rectangular grids (14/07/2018)
%
% Created on 25/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [P,it,res,flag] = SOR(P,RHS,U,V,GX,GY,imax,jmax,delx,dely,epsi,itermax,omg,nu,K2D,solver)

%------------------------- Iterative Solver -----------------------------%
	if solver == 1
		% Reshape P matrix into a column vector of stacked columns from rows of P, similar for RHS
		rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
        p = reshape(P(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);

        % Solve via Preconditioned Conjugate Gradient Algorithm without preconditioner
%         [p,flag,res,it] = pcg(K2D,rhs,epsi,itermax,[],[],p);

        % Alternatively use SOR-method (might be slower)
        [p, res, it] = sor_algorithm(K2D,rhs,p,omg,itermax,epsi);
		
		% Write pressure for next time step back into the matrix P
        P(2:imax,2:jmax) = reshape(p(:),[(imax-1) (jmax-1)]);
        
		clear p rhs; % Free memory
        
%--------------------------- Direct Solver ------------------------------%
    elseif solver == 2
		% Reshape P matrix into a column vector of stacked columns from rows of P, similar for RHS
		rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
        
        % Solve the linear system using a direct solver
        p = K2D\rhs;
        res = norm(K2D*p-rhs,2);
        it=1;   % Unused because direct solver
        
		% Write pressure for next time step back into the matrix P
        P(2:imax,2:jmax) = reshape(p(:),[(imax-1) (jmax-1)]);
        
		clear p rhs; % Free memory
        
%------------------------- Grid-oriented SOR ----------------------------%        
    else
        it = 1;
        res = 1; % Some value > epsi
        while( res > epsi && it < itermax )
			for i = 2:imax
				for j = 2:jmax
                    % Evaluate whether an interior, boundary or corner cell
                    % is computed
                    [c,Px,Py] = evaluate_P_terms(P,i,j,imax,jmax,delx,dely);
                    P(i,j) = (1-omg)*P(i,j) + omg*( c*( RHS(i,j) + Px + Py ) );
				end
            end
            p = reshape(P(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
            rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
	        res = norm(K2D*p-rhs,2);
			it = it+1;
        end
        clear p rhs; % Free memory
        
    end
    
end
	
	