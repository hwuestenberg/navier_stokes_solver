%######################################################
% Interpolate the velocity u in the real system from four positions on the 
% staggered grid
% Last: Created
%
% Created on 15/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [u,v] = interpolation(U,V,imax,jmax,delx,dely,x,y)
    % Ignore particles, that have been removed
    if( x <= 0 || y <= 0 )
        u = 0;
        v = 0;
        return;
    end
%------------------------ Horizontal Velocity ----------------------------%
    % Find cell, initial positions are allowed from x=[0.1,3.3] and y=[0,3.1]
    i = round(x/delx)+1;
    j = round((y+dely/2)/dely)+1;       % May change to symbolic integals using imax/jmax
    % Determine edge coordinates
    x1 = (i-1)*delx;
    x2 = i*delx;
    y1 = ((j-1)-0.5)*dely;
    y2 = (j-0.5)*dely;
    % Determine velocities at edges
    u1 = U(i-1,j-1);
    u2 = U(i,j-1);
    u3 = U(i-1,j);
    u4 = U(i,j);
    % Bilinearly interpolate u
    u = 1/delx/dely*( (x2-x)*(y2-y)*u1 + (x-x1)*(y2-y)*u2 +...
        (x2-x)*(y-y1)*u3 + (x-x1)*(y-y1)*u4 );
    
%------------------------- Vertical Velocity -----------------------------%
    % Find cell
    i = round((x+delx/2)/delx)+1;
    j = round(y/dely)+1;                % May change to symbolic integals using imax/jmax
    % Determine edge coordinates
    x1 = ((i-1)-0.5)*delx;
    x2 = (i-0.5)*delx;
    y1 = (j-1)*dely;
    y2 = j*dely;
    % Determine velocities at edges
    v1 = V(i-1,j-1);
    v2 = V(i,j-1);
    v3 = V(i-1,j);
    v4 = V(i,j);
    % Bilinearly interpolate u
    v = 1/delx/dely*( (x2-x)*(y2-y)*v1 + (x-x1)*(y2-y)*v2 +...
        (x2-x)*(y-y1)*v3 + (x-x1)*(y-y1)*v4 );
end

