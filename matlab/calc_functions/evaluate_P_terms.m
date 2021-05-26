%######################################################
% Evaluate presszre terms for the grid-oriented SOR method depending on 
% interior, boundary or corner cells
% Last: Created
%
% Created on 14/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [c,Px,Py] = evaluate_P_terms(P,i,j,imax,jmax,delx,dely)
    % Interior cell by default
    c = (delx^2*dely^2)/(2*delx^2+2*dely^2);
    Px = 1/delx^2*( P(i+1,j)+P(i-1,j) );
    Py = 1/dely^2*( P(i,j+1)+P(i,j-1) );
    
    % Left boundary
    if( i==2 )
        if( j==2 )          % Left and lower corner
            c = (delx^2*dely^2)/(delx^2+dely^2);  
            Px = 1/delx^2*( P(i+1,j) );
            Py = 1/dely^2*( P(i,j+1) );
        elseif( j==jmax )   % Left and upper corner
            c = (delx^2*dely^2)/(delx^2+dely^2); 
            Px = 1/delx^2*( P(i+1,j) );
            Py = 1/dely^2*( P(i,j-1) );
        else                % Left boundary
            c = (delx^2*dely^2)/(delx^2+2*dely^2);
            Px = 1/delx^2*( P(i+1,j) );
            Py = 1/dely^2*( P(i,j+1)+P(i,j-1) );
        end
    end
    
    % Right boundary
    if( i==imax )
        if( j==2 )          % Right and lower corner
            c = (delx^2*dely^2)/(delx^2+dely^2); 
            Px = 1/delx^2*( P(i-1,j) );
            Py = 1/dely^2*( P(i,j+1) );
        elseif( j==jmax )   % Right and upper corner
            c = (delx^2*dely^2)/(delx^2+dely^2);
            Px = 1/delx^2*( P(i-1,j) );
            Py = 1/dely^2*( P(i,j-1) );
        else                % Right boundary
            c = (delx^2*dely^2)/(delx^2+2*dely^2);
            Px = 1/delx^2*( P(i-1,j) );
            Py = 1/dely^2*( P(i,j+1)+P(i,j-1) );
        end
    end
    
    % Lower boundary
    if( j==2 )
        c = (delx^2*dely^2)/(2*delx^2+dely^2);
        Px = 1/delx^2*( P(i+1,j)+P(i-1,j) );
        Py = 1/dely^2*( P(i,j+1) );
    end
    
    % Upper boundary
    if( j==jmax )
        c = (delx^2*dely^2)/(2*delx^2+dely^2);
        Px = 1/delx^2*( P(i+1,j)+P(i-1,j) );
        Py = 1/dely^2*( P(i,j-1) );
    end
end



