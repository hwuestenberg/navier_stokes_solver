%######################################################
% Apply a Preconditioned Conjugate Gradient method to solve the pressure
% equation
%
%
% Created on 25/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [P,it,res,t_solver,resvec] = PCG(P,RHS,K2D,imax,jmax,epsi,itermax)



%------------------------- PCG Solver -----------------------------%
		
        % Reshape P matrix into a column vector of stacked columns from rows of P, similar for RHS
		rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
        p = reshape(P(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);

        
        % Solve via Preconditioned Conjugate Gradient Algorithm without preconditioner
        t_solver = tic;
        
        [p,flag,res,it,resvec] = pcg(K2D,rhs,epsi,itermax,[],[],p);
		
        t_solver = toc(t_solver);
        
%         % Plot residuum norm over all iterations
%         semilogy(1:it,resvec(1:it));
%         ylim([1e-3 1e+1])
%         xlim([0 20])
%         drawnow;
        
        it = it+1;
		% Write pressure for next time step back into the matrix P
        P(2:imax,2:jmax) = reshape(p(:),[(imax-1) (jmax-1)]);
        
		clear p rhs; % Free memory
    
    
end
	
	