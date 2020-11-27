%######################################################
% Visualize data of the flow simulation using a mesh-plot
% How to use 'vis_mode' parameter:
% '-1' Off
% '0' Multiple figures
% '1' subplots
% '2' Velocity field & Pressure & Streamlines
% '3' Velocity field & Vorticity & Streamlines
%
% Created on 27/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function visual(U,V,P,Psi,Zeta,res,vis_mode,imax,jmax,problem)
if vis_mode == -1
    return;
end
%--------------------------- Transpose Data -----------------------------%
% In order to visualise data correctly, one needs to transose each matrix
% U=U'; V=V'; P=P'; Psi=Psi'; Zeta=Zeta';

%---------------------------- Axis Limits -------------------------------%

    % Estimate limits for each physical quantity (does not work yet)
%     c = 0.5;    % Correction for limits
%     maxU = round(max(max(U))*(1+c),1); minU = round(min(min(U))*(1-c),1);
%     maxV = round(max(max(V))*(1+c),1); minV = round(min(min(V))*(1-c),1);
%     maxP = round(max(max(P))*(1+c),1); minP = round(min(min(P))*(1-c),1);
%     maxPsi = round(max(max(Psi))*(1+c),1); minPsi = round(min(min(Psi))*(1-c),1);
%     maxZeta = round(max(max(Zeta))*(1+c),1); minZeta = round(min(min(Zeta))*(1-c),1);
%     % Check limits for double zeros
%     if maxU < 1e-10; maxU = 1e-10; end
%     if maxV < 1e-10; maxV = 1e-10; end
%     if maxP < 1e-10; maxP = 1e-10; end
%     if maxPsi < 1e-10; maxPsi = 1e-10; end
%     if maxZeta < 1e-10; maxZeta = 1e-10; end

    % Set experience-based values for axis limits
    switch(problem)
        case 1
            maxU = 2; minU = -maxU;
            maxV = 1; minV = -maxV;
            maxP = 5; minP = -maxP;
            maxPsi = 0.3; minPsi = -maxPsi;
            maxZeta = 4; minZeta = -maxZeta;
        case 2
            maxU = 2.5; minU = -maxU;
            maxV = 2.5; minV = -maxV;
            maxP = 3; minP = -maxP;
            maxPsi = 10; minPsi = -maxPsi;
            maxZeta = 6; minZeta = -maxZeta;
        case 3
            maxU = 2; minU = -maxU;
            maxV = 1; minV = -maxV;
            maxP = 5; minP = -maxP;
            maxPsi = 0.3; minPsi = -maxPsi;
            maxZeta = 4; minZeta = -maxZeta;
        case 4
            maxU = 2; minU = -maxU;
            maxV = 1; minV = -maxV;
            maxP = 5; minP = -maxP;
            maxPsi = 0.3; minPsi = -maxPsi;
            maxZeta = 4; minZeta = -maxZeta;    
    end
    
%------------------------------- Plots ---------------------------------    
    % Multiple plots
    switch vis_mode
        case 0
            figU = figure(1);
            mesh(U(2:imax-1,2:jmax)')
            title('Velocity U')
            xlabel('x')
            ylabel('y')
            zlabel('U [m/s]')
            xlim([1 imax+1])
            ylim([1 jmax+1])
            zlim([minU maxU])

            figV = figure(2);
            mesh(V(2:imax,2:jmax-1)')
            title('Velocity V')
            xlabel('x')
            ylabel('y')
            zlabel('V [m/s]')
            xlim([1 imax+1])
            ylim([1 jmax+1])
            zlim([minV maxV])

            figP = figure(3);
            mesh(P(2:imax,2:jmax)')
            title('Pressure P')
            xlabel('x')
            ylabel('y')
            zlabel('P [N/m^2]')
            xlim([1 imax+1])
            ylim([1 jmax+1])
            zlim([minP maxP])
            

        % Single subplot
        case 1
            figure(1)
            figU = subplot(2,3,1);
            mesh(U(2:imax-1,2:jmax)')
            title('Velocity U')
            xlabel('x')
            ylabel('y')
            zlabel('U [m/s]')
            xlim([1 imax+1])
            ylim([1 jmax+1])
            zlim([minU maxU])

            figV = subplot(2,3,2);
            mesh(V(2:imax,2:jmax-1)')
            title('Velocity V')
            xlabel('x')
            ylabel('y')
            zlabel('V [m/s]')
            xlim([1 imax+1])
            ylim([1 jmax+1])
            zlim([minV maxV])

            figP = subplot(2,3,3);
            mesh(P(2:imax,2:jmax)')
            title('Pressure P')
            xlabel('x')
            ylabel('y')
            zlabel('P [N/m^2]')
            xlim([1 imax+1])
            ylim([1 jmax+1])
            zlim([minP maxP])

            figPsi = subplot(2,3,4);      
            contour(Psi')
            hold on
            quiver(2:imax,2:jmax,U(2:imax,2:jmax)',V(2:imax,2:jmax)')
            title('Streamlines')
            xlabel('x')
            ylabel('y')
            hold off
            
            figRes = subplot(2,3,6);
            semilogy(1:length(res),res,'--or')
            title('Convergence')
            xlabel('Time step')
            ylabel('Residuum')
         
    % Velocity field with quiver, Pressure levels with contour and
    % Streamlines
        case 2
            figPsi = subplot(1,2,1);
            contour(Psi')
            hold on
            quiver(2:imax,2:jmax,U(2:imax,2:jmax)',V(2:imax,2:jmax)')
            title('Streamlines')
            xlabel('x')
            ylabel('y')            
            hold off
            
            figP = subplot(1,2,2);
            contour(P')
            title('Pressure levels')
            xlabel('x')
            ylabel('y')
            hold on
            quiver(2:imax,2:jmax,U(2:imax,2:jmax)',V(2:imax,2:jmax)')
            colorbar
            hold off
    
    % Velocity field with quiver, Streamlines and Vorticity with contour
        case 3
            figPsi = subplot(1,2,1);      
            contour(Psi')
            hold on
            quiver(2:imax,2:jmax,U(2:imax,2:jmax)',V(2:imax,2:jmax)')
            title('Streamlines')
            xlabel('x')
            ylabel('y')            
            hold off

            figZeta = subplot(1,2,2);
            contour(Zeta')
            hold on
            quiver(2:imax,2:jmax,U(2:imax,2:jmax)',V(2:imax,2:jmax)')
            title('Vorticity')
            xlabel('x')
            ylabel('y')
            hold off
            
    end
    
    % Save picture of current result
    if problem == 1
        problem_str = 'DrivenCavity_';
    elseif problem == 2
        problem_str = 'KarmanVortexStreet_';
    elseif problem == 3
        problem_str = 'FlowAboveStair_';
    end
    str = strcat('results/',problem_str,num2str(imax-1),'x',num2str(jmax-1),'_',num2str(round(1000*rand(1),0)));
    print(str,'-dpng');
end