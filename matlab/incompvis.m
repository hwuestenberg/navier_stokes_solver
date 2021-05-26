%######################################################
% This script executes a Navier-Stokes-based simulation for incompressible,
% viscid, steady and two-dimensional flow. It is required to run the script
% progpar in before to intialize a problem together with its parameters and
% a solver.
%
% Great changes:
% Added individual boundary configuration Free-/No-slip (28/06/2018)
% Added Multigrid solver for V/W-cycle and Full Multigrid
%
% Current Issues: 
% Grid-wise solver has an offset error in residuum (14/07/2018)
%
% Created on 27/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
%% Initilization
load('sim_data_init'); % Load data set

% Change default figure size and position
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'DefaultAxesFontSize',20)
set(groot, 'defaultFigurePosition',[0 0 1 1])

%% Run Simulation
% Check stability condition 26.1 and, eventually, decrease time step size
stab1(delt,delx,dely,nu);


[U,V,P] = initgrid(imax,jmax,U_I,V_I,P_I);  % Set initial values for U,V,P
[K2D] = init_linSys(imax,jmax,delx,dely);   % Set linear system for SOR



% Visualisation of Streaklines and Particle tracers
if visualise
    part_x_streak = zeros(N,T_end/(delt_n*delt)); % Vectors to store injected particles after 4 time steps
    part_y_streak = zeros(N,T_end/(delt_n*delt)); % Coordinates for injected particle tracing
    part_x_trace = zeros(N,T_end/delt);     % Vectors to store particle ...
    part_y_trace = zeros(N,T_end/delt);     % Coordinates for tracing
    [part_x_trace(:,1),part_y_trace(:,1)] = set_particles(N,ug,og,delx);
    last_part_x = part_x_trace(:,1);
    last_part_y = part_y_trace(:,1);
end


% Recording of individual frames for an animation
if movie_mode
    M = struct('cdata',[],'colormap',[]);   % Auxiliary struct for frames
end


% Struct containing:
% Residuals
% Time measurement for solver
% Total simulated time
% Iteration counter
% Vector containing residuals of each iteration
solv = struct('res',zeros(T_end/delt,1), 'time',zeros(T_end/delt,1), 'T',0, 'iter',zeros(T_end/delt,1), 'resvec',zeros(T_end/delt,1));                         

i = 1;  % Iteration and Frame counter
t_loop = tic;
while solv.T < (T_end-delt/2)        % Avoid failure due to machine accuracy
    % Set boundary conditions
    [U,V,P,F,G] = boundary_values(U,V,P,imax,jmax,wl,wr,wt,wb);
    [U,V,P] = spec_boundary_values(U,V,P,imax,jmax,problem);
    
    
    % Compute right-hand side 
    [F,G] = calc_FG(U,V,F,G,imax,jmax,delt,delx,dely,GX,GY,alpha,nu);
    RHS = calc_RHS(F,G,imax,jmax,delt,delx,dely);
    
    
    
    % Solve the pressure equation
    
    % Grid-oriented SOR
    if solver == 0
        [P,solv.iter(i),solv.res(i),solv.time(i),solv.resvec] = Grid_SOR(P,RHS,K2D,imax,jmax,epsi,omg,itermax,delx,dely);
        
    % SOR
    elseif solver == 1
        [P,solv.iter(i),solv.res(i),solv.time(i),solv.resvec] = SOR(P,RHS,K2D,imax,jmax,epsi,omg,itermax);
    
    % PCG
    elseif solver == 2
        [P,solv.iter(i),solv.res(i),solv.time(i),solv.resvec] = PCG(P,RHS,K2D,imax,jmax,epsi,itermax);
    
    % Multigrid
    elseif solver > 2 && solver < 7
        [P,solv.iter(i),solv.res(i),solv.time(i)] = multigrid(P,RHS,K2D,imax,jmax,delx,dely,epsi,omg,solver);
    
    % Direct Solver
    elseif solver == 7
        [P,solv.iter(i),solv.res(i),solv.time(i)] = Direct_Solver(P,RHS,K2D,imax,jmax);
        
    end
    
    
    
    
    
    % Compute new velocities
    [U,V] = calc_UV(U,V,P,F,G,imax,jmax,delt,delx,dely);
    
    
    % Check stability and, if necessary, stop simulation
    stabCond = stab23(U,V,delt,delx,dely,solv.T);   % Check CFL condition 26.2/26.3
    if(~stabCond)
        break;
    end
    
    
    % Compute Streaklines and Particle tracers
    if visualise
        [part_x_trace(:,i),part_y_trace(:,i)] = particletrace(U,V,imax,jmax,delx,dely,delt,last_part_x,last_part_y,problem);
        [part_x_streak,part_y_streak] = streaklines(U,V,imax,jmax,delx,dely,delt,part_x_streak,part_y_streak,N,ug,og,delt_n,solv.T,T_end,problem);
    end
    
    
    % Time stepping
    solv.T = solv.T + delt;
    
    
    % Save each frame for an animation
    if movie_mode
        [Psi,Zeta] = psi_zeta(U,V,imax,jmax,delx,dely);
        visual(U,V,P,Psi,Zeta,solv.res,vis_mode,imax,jmax,problem)
        drawnow
        M(i) = getframe(gcf);
    end
    
    
    i = i+1;    % Increase Iteration and Frame counter
    
    
    if visualise
        % Save last coordinates for particle trace function
        last_part_x = part_x_trace(:,i-1);
        last_part_y = part_y_trace(:,i-1);
    end
end
t_loop = toc(t_loop);


% Compute Vorticity and Stream function
[Psi,Zeta] = psi_zeta(U,V,imax,jmax,delx,dely);


% Visualise the results
visual(U,V,P,Psi,Zeta,solv.res,vis_mode,imax,jmax,problem);


% Visualise Streaklines and Particle tracers
if visualise
    particlevis(part_x_trace,part_y_trace,part_x_streak,part_y_streak,delx,dely,imax,jmax,problem);
end


% Visualise the results as movie
if movie_mode
    movie(gcf,M,2)
end


% Time measurment and Covergence behaviour
if solver_measurement
    evaluate_Solver(solv,epsi,omg,solver);
end


fprintf('The simulation finished after %4.3f seconds\n',t_loop);


% Save results
save('sim_data_result')


