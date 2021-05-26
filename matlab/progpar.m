%######################################################
% Configure parameters for the Navier-Stokes solver
%
% Created on 24/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
clear all;

% Choose method and parameters
solver = 6;     % '0' Grid-oriented | '1' SOR | '2' PCG
                % '3' Simple Multigrid | '4' V-Cycle Multigrid 
                % '5' W-Cycle Multigrid | '6' Full Multigrid
                % '7' Direct
                
vis_mode = 1;               % See function description in 'visual.m'
visualise = 0;              % Activate Streaklines and Particle tracing
solver_measurement = 0;     % Activate plotting of solver information
movie_mode = 0;             % Activate recording

% Choose the problem (Check for corresponding boundary conditions)
% 1 = Driven Cavity
% 2 = Van Karman vortex street
% 3 = Flow above a Stair
% 4 = Tunnel with Constricitons
problem = 2;

switch(problem)
%----------------------------- Driven Cavity -----------------------------
    case 1
        % Discretisation
        imax = 16;
        jmax = 16;
        imax = imax+1;  % Correction because Matlab-Arrays start with 1 not 0
        jmax = jmax+1;
        delx = 0.2;
        dely = 0.2;
        delt = 0.02;
        T_end = 3;

        % Iterative solving
        itermax = 150;
        epsi = 0.01;
        omg = 1.7;
        alpha = 0.12;

        % Fluid properties, initial velocity and pressure field
        nu = 0.4;
        GX = 0;
        GY = 0;
        U_I = 0;
        V_I = 0;
        P_I = 0;

        % Geometry definition
        % 1 = Free slip
        % 2 = No-slip
        % 3 = Inflow/Outflow (Inflow at left wall)
        wl = 2;     % Left wall
        wr = 2;     % Right wall
        wt = 2;     % Top wall
        wb = 2;     % Bottom wall
        
        % Streaklines and Particle Tracing
        N = 16;
        ug = 2*dely;
        og = 14*dely;
        delt_n = 4;
        
%------------------------ Van Karman Vortex Street -----------------------
    case 2
        % Discretisation
        imax = 64;
        jmax = 16;
        imax = imax+1;  % Correction because Matlab-Arrays start with 1 not 0
        jmax = jmax+1;
        delx = 0.2;
        dely = 0.2;
        delt = 0.02;
        T_end = 20;

        % Iterative solving
        itermax = 150;
        epsi = 0.01;
        omg = 1.7;
        alpha = 0.9;

        % Fluid properties, initial velocity and pressure field
        nu = 1.3414e-5;
        GX = 0;
        GY = 0;
        U_I = 1;
        V_I = 0;
        P_I = 0;

        % Geometry definition
        % 1 = Free slip
        % 2 = No-slip
        % 3 = Inflow/Outflow (Inflow at left wall)
        wl = 3;     % Left wall
        wr = 3;     % Right wall
        wt = 1;     % Top wall
        wb = 1;     % Bottom wall        
        
        % Streaklines and Particle Tracing
        N = 16;
        ug = 2*dely;
        og = 14*dely;
        delt_n = 4;

%--------------------------- Flow above a Stair --------------------------
    case 3
        % Discretisation
        imax = 64;
        jmax = 16;
        imax = imax+1;  % Correction because Matlab-Arrays start with 1 not 0
        jmax = jmax+1;
        delx = 0.2;
        dely = 0.2;
        delt = 0.02;
        T_end = 20;

        % Iterative solving
        itermax = 150;
        epsi = 0.01;
        omg = 1.7;
        alpha = 0.12;

        % Fluid properties, initial velocity and pressure field
        nu = 0.1;
        GX = 0;
        GY = 0;
        U_I = 1;
        V_I = 0;
        P_I = 0;

        % Geometry definition
        % 1 = Free slip
        % 2 = No-slip
        % 3 = Inflow/Outflow (Inflow at left wall)
        wl = 3;     % Left wall
        wr = 3;     % Right wall
        wt = 2;     % Top wall
        wb = 2;     % Bottom wall
        
        % Streaklines and Particle Tracing
        N = 16;
        ug = 2*dely;
        og = 14*dely;
        delt_n = 4;
        
%------------------------ Tunnel with Constricitons ----------------------
    case 4
        % Discretisation
        imax = 64;
        jmax = 16;
        imax = imax+1;  % Correction because Matlab-Arrays start with 1 not 0
        jmax = jmax+1;
        delx = 0.2;
        dely = 0.2;
        delt = 0.02;
        T_end = 0.5;

        % Iterative solving
        itermax = 150;
        epsi = 0.01;
        omg = 1.7;
        alpha = 0.9;

        % Fluid properties, initial velocity and pressure field
        nu = 1.3414e-5;
        GX = 0;
        GY = 0;
        U_I = 1;
        V_I = 0;
        P_I = 0;

        % Geometry definition
        % 1 = Free slip
        % 2 = No-slip
        % 3 = Inflow/Outflow (Inflow at left wall)
        wl = 3;     % Left wall
        wr = 3;     % Right wall
        wt = 1;     % Top wall
        wb = 1;     % Bottom wall
        
        % Streaklines and Particle Tracing
        N = 16;
        ug = 2*dely;
        og = 14*dely;
        delt_n = 4;

end

% Output current configuration
switch(solver)
    case 0
        solver_str = 'Grid-Oriented';
    case 1
        solver_str = 'SOR';
    case 2
        solver_str = 'PCG';
    case 3
        solver_str = 'Simple Multigrid';
    case 4
        solver_str = 'V-Cycle Multigrid';
    case 5
        solver_str = 'W-Cycle Multigrid';
    case 6
        solver_str = 'Full Multigrid';
    case 7
        solver_str = 'Direct';
end

switch(movie_mode)
    case 0
        movie_str = 'no';
    case 1
        movie_str = 'with';
end

switch(problem)
    case 0
        problem_str = 'Debug';
    case 1
        problem_str = 'Driven Cavity';
    case 2
        problem_str = 'Von Karman Vortex Street';
    case 3
        problem_str = 'Flow above a Stair';
    case 4
        problem_str = 'Tunnel with Constrictions';
end
fprintf('\nInitialized Parameters for %s problem with\n%s solver, visualization mode %i and %s movie\n',problem_str,solver_str,vis_mode,movie_str);
fprintf('Grid size is set to [%ix%i] (x by y)\n',imax-1,jmax-1);

% Save configuration for simulation (Overwriting old configuration)
clear problem_str movie_str solver_str
save('sim_data_init');

