%######################################################
% Set boundary conditions for a given problem using the 'problem' parameter
% Last: Added problem 4: Tunnel with constrictions (30/06/2018)
%
% Created on 28/06/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [U,V,P] = spec_boundary_values(U,V,P,imax,jmax,problem)
%-------------------------- Debugging Problem ---------------------------%
    if problem == 0
        U(1,2:jmax) = 2;    % const. horizontal velocity at left boundary
    end
%-------------------------- Driven Cavity -------------------------------%    
    if problem == 1
        U(2:imax,jmax+1) = 2-U(2:imax,jmax);    % Top with U=1 (avg.)
    end
%-------------------- Von Karmann Vortex Street -------------------------%
    if problem == 2
        % Inflow with linear velocity profile
        for j = 2:jmax
            U(1,j) = (j-1)/(jmax-1);
        end
        V(1,2:jmax) = -V(2,2:jmax); % V=0 at left boundary
        
        % Stumbling block in central and forward position of the grid
        yfifth = jmax/5;
        offset_x = 4*yfifth;
        
        left = round(offset_x)+1;
        right = round(offset_x+yfifth)+1;
        lower = round(2*yfifth)+1;
        upper = round(3*yfifth)+1;
        
        for i = left:right
            for j = lower:upper
                U(i,j) = 0;
                V(i,j) = 0;
                P(i,j) = 0;
            end
        end
        
        % Pressure conditions for internal walls
        % No pressure gradient across internal walls
        P(left,lower:upper) = P(left-1,lower:upper);    % Left
        P(right,lower:upper) = P(right+1,lower:upper);  % Right
        P(left:right,upper) = P(left:right,upper+1);    % Top
        P(left:right,lower) = P(left:right,lower-1);    % Bottom
    end
%------------------------ Flow above a Stair ----------------------------%
    if problem == 3
        yhalf = round((jmax-1)/2);
        % Inflow constant profile
        for j = yhalf+1:jmax
            U(1,j) = 1;
            V(1,j) = -V(2,j);
        end
        
        % Stair
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



