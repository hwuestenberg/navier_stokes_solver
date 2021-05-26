%######################################################
% Evaluate F and G terms for the grid-oriented SOR method depending on 
% interior, boundary or corner cells
% Last: Causes divergence for T>0.3 beginning at corners, why?! Therefore,
% it is not used at the moment.
%
% Created on 14/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [f,g] = evaluate_FG_terms(F,G,i,j,imax,jmax)
    % Interior cell by default
    f = F(i-1,j)-F(i,j);
    g = G(i,j-1)-G(i,j);
    
    % Left boundary
    if( i==2 )
        if( j==2 )          % Left and lower corner
            f = -F(i,j);
            g = -G(i,j);
        elseif( j==jmax )   % Left and upper corner
            f = -F(i,j);
            g = G(i,j-1);
        else                % Left boundary
            f = -F(i,j);
            g = G(i,j-1)-G(i,j);
        end
    end
    
    % Right boundary
    if( i==imax )
        if( j==2 )          % Right and lower corner
            f = F(i-1,j);
            g = -G(i,j);
        elseif( j==jmax )   % Right and upper corner
            f = F(i-1,j);
            g = G(i,j-1);
        else                % Right boundary
            f = F(i-1,j);
            g = G(i,j-1)-G(i,j);
        end
    end
    
    % Lower boundary
    if( j==2 )
        f = F(i-1,j)-F(i,j);
        g = -G(i,j);
    end
    
    % Upper boundary
    if( j==2 )
        f = F(i-1,j)-F(i,j);
        g = G(i,j-1);
    end
end



