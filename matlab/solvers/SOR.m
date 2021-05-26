%######################################################
% Apply a Successive over-relaxation (SOR) method to solve the pressure
% equation
%
%
% Created on 25/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [P,it,res,t_solver,resvec] = SOR(P,RHS,K2D,imax,jmax,epsi,omg,itermax)



%------------------------- SOR Solver -----------------------------%
		
        % Reshape P matrix into a column vector of stacked columns from rows of P, similar for RHS
		rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
        p = reshape(P(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);


        % Alternatively use SOR-method (might be slower)
        t_solver = tic;
        
        [p,it,res,resvec] = sor_algorithm(K2D,rhs,p,omg,itermax,epsi);
		
        t_solver = toc(t_solver);
        
        % Plot residuum norm over all iteraions
%         semilogy(1:it,resvec(1:it));
%         ylim([1e-3 1e+1])
%         xlim([0 100])
%         drawnow;
        
        
		% Write pressure for next time step back into the matrix P
        P(2:imax,2:jmax) = reshape(p(:),[(imax-1) (jmax-1)]);
        
		clear p rhs; % Free memory
    
    
end
	
	