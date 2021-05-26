%######################################################
% Apply a variety of multigrid methods to efficiently solve the pressure 
% equation
%
%
% Last: Created
%
%
% Created on 10/08/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [P,it,res,t_solver] = multigrid(P,RHS,K2D,imax,jmax,delx,dely,epsi,omg,solver)


    % Reshape P matrix into a column vector of stacked columns from rows of P, similar for RHS
    rhs = reshape(RHS(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
    p = reshape(P(2:imax,2:jmax),[(imax-1)*(jmax-1) 1]);
    
    
    L = 0;          % Initialize grid level to zero
    
    iter1 = 20;     % Set smoothing iterations
    iter2 = iter1;

    
    
%------------------- Multigrid Solver without recursion ------------------%
     if solver == 3
        
        % Determine size of x vector for local level
		loc_size = (imax-1)*(jmax-1)/4^(L+1);
		loc_i = (imax-1)/2^(L+1);	% local grid size in x-direction
		loc_j = (jmax-1)/2^(L+1);	% local grid size in y-direction
        
        p = sor_algorithm(K2D,rhs,p,omg,iter1,epsi);  % Smoothing applied iter1 times
		
		
        res_fine = rhs-K2D*p;                                					% Compute residual
        res_coarse = restriction(res_fine,imax,jmax,L);        					% Restriciton
        K2D_coarse = init_linSys(loc_i+1,loc_j+1,delx*2^(L+1),dely*2^(L+1));	% Compute coarse matrix
		
        p_coarse = K2D_coarse\res_coarse;			% Direct solver for coarse grid
        
        p = p + prolongation(p_coarse,imax,L);      % Prolongation
        
		
		[p,res,it]  = sor_algorithm(K2D,rhs,p,omg,iter2,epsi);  % Smoothing applied iter2 times
        
     end
    
     
     
     
%----------------------- V-Cycle Multigrid Solver ------------------------%
     if solver == 4
         t_solver = tic;
         [p,it,res] = VCycleMultigrid(p,rhs,K2D,imax,jmax,delx,dely,epsi,omg,iter1,iter2,L);
         t_solver = toc(t_solver);
     end
     
     
     
     
%----------------------- W-Cycle Multigrid Solver ------------------------%
     if solver == 5
        t_solver = tic;
        [p,it,res] = WCycleMultigrid(p,rhs,K2D,imax,jmax,delx,dely,epsi,omg,iter1,iter2,L);
        t_solver = toc(t_solver);
        
     end
    
     
     
%------------------------- Full Multigrid Solver -------------------------%
     if solver == 6
        t_solver = tic;
        [p,it,res] = FullMultigrid(p,rhs,K2D,imax,jmax,delx,dely,epsi,omg,iter1,iter2,L);
        t_solver = toc(t_solver);
        
     end
     
     
     
    P(2:imax,2:jmax) = reshape(p(:),[(imax-1) (jmax-1)]);   % Backtransform into a matrix
    clear p rhs % Free memory
     
end
	
	