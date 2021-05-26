%######################################################
% Apply a grid-oriented Successive over-relaxation (SOR) method to solve 
% the pressure equation
%
%
% Created on 25/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [P,it,res,t_solver,resvec] = Grid_SOR(P,RHS,K2D,imax,jmax,epsi,omg,itermax,delx,dely)



%------------------------- Grid-oriented SOR ----------------------------%        
    it = 1;
    res = 1; % Some value > epsi
    
    resvec = zeros(itermax,1);
    
    while( res > epsi && it < itermax )
        t_solver = tic;
        
        for i = 2:imax
            for j = 2:jmax
                % Evaluate whether an interior, boundary or corner cell
                % is computed
                [c,Px,Py] = evaluate_P_terms(P,i,j,imax,jmax,delx,dely);
                P(i,j) = (1-omg)*P(i,j) + omg*( c*( RHS(i,j) + Px + Py ) );
            end
        end
        
        t_solver = toc(t_solver);
        
        p = reshape(P(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
        rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);

        res = norm(K2D*p-rhs,2);
        resvec(it) = res;
        it = it+1;
    end
    
    % Free memory
    clear p rhs; 
    
    
end
	
	