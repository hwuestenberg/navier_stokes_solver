%######################################################
% A W-cycle Multigrid solver to effciently solve linear systems
%
%
% Last: Created
%
%
% Created on 13/08/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [x,it,res] = WCycleMultigrid(x,b,A,size_row_p1,size_col_p1,delx,dely,epsi,omg,iter1,iter2,L)
     
%----------------------- W-Cycle Multigrid Solver ------------------------%
    % Determine size of x vector for local level
    loc_size = (size_row_p1-1)*(size_col_p1-1)/4^(L+1);
	loc_i = (size_row_p1-1)/2^(L+1);	% local grid size in x-direction
	loc_j = (size_col_p1-1)/2^(L+1);	% local grid size in y-direction

    x = sor_algorithm(A,b,x,omg,iter1,epsi);    % 1st Smoothing applied iter1 times


    res_coarse = b-A*x;   												% Compute residual
    res_coarse = restriction(res_coarse,size_row_p1,size_col_p1,L);		% Restriction
	A_coarse = init_linSys(loc_i+1,loc_j+1,delx*2^(L+1),dely*2^(L+1));	% Compute coarse matrix
	
    % Check level, apply direct solver on coarsest level
    if L == 0
        x_coarse = WCycleMultigrid(zeros(loc_size,1),res_coarse,A_coarse,size_row_p1,size_col_p1,delx,dely,epsi,omg,iter1,iter2,L+1);
    elseif L > 0 && L < 2
        % Recursive call with x=0, b=res_coarse, A=A_coarse L=L+1
        x_coarse = WCycleMultigrid(zeros(loc_size,1),res_coarse,A_coarse,size_row_p1,size_col_p1,delx,dely,epsi,omg,iter1,iter2,L+1);
        % Recursive call with x=x_coarse, b=res_coarse, A=A_coarse L=L+1
		x_coarse = WCycleMultigrid(x_coarse,res_coarse,A_coarse,size_row_p1,size_col_p1,delx,dely,epsi,omg,iter1,iter2,L+1);

    else
        % Solve directly on coarsest grid
        x_coarse = A_coarse\res_coarse;

    end
    
    
    % Correct x vector
    x = x + prolongation(x_coarse,size_row_p1,L);


    [x,it,res] = sor_algorithm(A,b,x,omg,iter2,epsi);  % Smoothing applied iter2 times
        
end