%######################################################
% Call the particletrace subroutine, that propagates all particles and
% detects when an particle leaves the grid
%
% Last: Set removing conditions for problem 1&2 (16/07/2018)
%
% Created on 16/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [part_x,part_y] = particletrace(U,V,imax,jmax,delx,dely,delt,part_x,part_y,problem)
    % Compute new position for next time_step
    [part_x,part_y] = new_position(U,V,imax,jmax,delx,dely,delt,part_x,part_y);
    
    % Detect particles leaving the grid
    for i = 1:length(part_x)
        % x-direction
        if( part_x(i) < 0.5*delx || part_x(i) > (imax-0.5)*delx )     % Eventually decrease delx to delx/2 and imax-0.5
            part_x(i) = -1; % Set as removed
            part_y(i) = -1;
        % y-direction    
        elseif( part_y(i) < 0.5*dely || part_y(i) > (jmax-0.5)*dely ) % Same as above
            part_x(i) = -1; % Set as removed
            part_y(i) = -1;
        end
    end
    
    % Detect particles that contact internal walls (problem-specific)
%-------------------------- Driven Cavity -------------------------------%    
    if problem == 1
        % No internal walls
    end
%-------------------- Von Karmann Vortex Street -------------------------%
    if problem == 2
        % Stumbling block in central and forward position of the grid
        yfifth = jmax/5;
        offset_x = 4*yfifth;
        left = round(offset_x)+1;
        right = round(offset_x+yfifth)+1;
        lower = round(2*yfifth)+1;
        upper = round(3*yfifth)+1;
        
        for i = 1:length(part_x)
            if( (part_x(i) > (left-1)*delx && part_x(i) < (right+1)*delx) && ...
                   part_y(i) > (lower-1)*dely && part_y(i) < (upper+1)*dely )       
               % Eventually decrease delx to delx/2 and imax-0.5 
               % Also, correct interpolation function
                part_x(i) = -1; % Set as removed
                part_y(i) = -1;
            end
        end
    end
%------------------------ Flow above a Stair ----------------------------%
    if problem == 3
        yhalf = round((jmax-1)/2);
        
        % Geometry within the flow field
        % Horizontal velocity
        U(1:jmax,2:yhalf) = 0;
        % Vertical velocity
        V(2:jmax,1:yhalf) = 0;
    end
%---------------------- Tunnel with Constrictions ------------------------%
    if problem == 4
        yfifth = round((jmax-1)/5); % Scaling that fits the given figure
        xthird = round((imax-1)/3);
        
        % Inflow constant profile
        U(1,yfifth*2+1:yfifth*3+1) = 1;
        V(1,yfifth*2+1:yfifth*3+1) = -V(2,yfifth*2+1:yfifth*3+1);
        
        % Geometry within the flow field
        % Forward constriction
        U(1:xthird+1,1:2*yfifth+1) = 0;
        U(1:xthird+1,3*yfifth+1:jmax) = 0;
        % Rear constriction
        U(2*xthird+1:3*xthird+1,1:2*yfifth+1) = 0;
        U(2*xthird+1:3*xthird+1,3*yfifth+1:jmax) = 0;
        
        % Forward constriction
        V(1:xthird+1,1:2*yfifth+1) = 0;
        V(1:xthird+1,3*yfifth:jmax) = 0;
        % Rear constriction
        V(2*xthird+1:3*xthird+1,1:2*yfifth+1) = 0;
        V(2*xthird+1:3*xthird+1,3*yfifth:jmax) = 0;
    end    
end

