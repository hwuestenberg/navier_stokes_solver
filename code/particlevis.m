%######################################################
% Plot Particle Tracer positions at every time step and Streaklines
% Last: Set removing conditions for problem 1&2 (16/07/2018)
%
% Created on 16/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function particlevis(part_x_trace,part_y_trace,delx,dely,imax,jmax,problem)
    % Setup an vector of different colours
    b = 'ob';
    r = 'or';
    % Create figure
    figure(10)
    hold on
    % Plot dots for each time step
    for i = 1:size(part_x_trace,2)
        if mod(i,2)
            plot(part_x_trace(:,i),part_y_trace(:,i),b);
        else
            plot(part_x_trace(:,i),part_y_trace(:,i),r);
        end
    end
    % Plot Boundaries
%-------------------------- Grid Boundaries ------------------------------%
    plot([0 imax*delx],[0 0],'-k',[0 imax*delx],[jmax*dely jmax*dely],'-k',[0 0],[0 jmax*dely],'-k',[imax*delx imax*delx],[0 jmax*dely],'-k')
%--------------------------- Driven Cavity -------------------------------%
    if problem == 1
        % No internal walls
    end
%--------------------- Von Karmann Vortex Street -------------------------%
    if problem == 2
        yfifth = jmax/5;
        offset_x = 4*yfifth;
        left = (round(offset_x)+1)*delx;
        right = (round(offset_x+yfifth)+1)*delx;
        lower = (round(2*yfifth)+1)*dely;
        upper = (round(3*yfifth)+1)*dely;
        block = polyshape([left left right right],[upper lower lower upper]);
        plot(block);
    end
    % Set axis parameters
    title('Particle Tracer');
    xlabel('x [m]')
    ylabel('y [m]')
    xlim([0 imax*delx])
    ylim([0 imax*delx])
end