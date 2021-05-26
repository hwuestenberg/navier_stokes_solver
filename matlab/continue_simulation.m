%######################################################
% Continues a fluid simulation that was run using the incompvis script.
% Parameters for the simulation are loaded using the progpar script.
% Last: Renamed choose variable 'fast' to 'solver'
%
% Created on 30/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
%% Initilization
load('sim_data_result'); % Load data set

% Change default figure size and position
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0.25 1 0.5])

%% Run Simulation
% Check stability condition 26.1 and, eventually, decrease time step size
stab1(delt,delx,dely,nu);

res = zeros(10,1);  % Save residuals
i=1;                % Iteration counter
T_end   = T_end+10*delt;
tic;
while T < (T_end-delt/2)
    [U,V,P,F,G] = boundary_values(U,V,P,imax,jmax,wl,wr,wt,wb);
    [U,V,P] = spec_boundary_values(U,V,P,imax,jmax,problem);
    
    [F,G] = calc_FG(U,V,F,G,imax,jmax,delt,delx,dely,GX,GY,alpha,nu);
    RHS = calc_RHS(F,G,imax,jmax,delt,delx,dely);
    
    [P,it,res(i)] = SOR(P,RHS,U,V,GX,GY,imax,jmax,delx,dely,epsi,itermax,omg,nu,K2D,solver);
    
    [U,V] = calc_UV(U,V,P,F,G,imax,jmax,delt,delx,dely);
    stab23(U,V,delt,delx,dely,T);   % Check CFL conditions 26.2 and 26.3 
    T = T+delt;                     % Time stepping
    i = i+1;
end
t_loop = toc;

% Compute Vorticity and Stream function
[Psi,Zeta] = psi_zeta(U,V,imax,jmax,delx,dely);

% Visualize the results
visual(U,V,P,Psi,Zeta,res,vis_mode,imax,jmax,problem);

fprintf('The simulation finished after %4.3f seconds\n',t_loop);
% Save results
save('sim_data_result')
