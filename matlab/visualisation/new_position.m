%######################################################
% Compute new positions for each particle based on there current velocity
% vector
% Last: Created
%
% Created on 15/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [part_x,part_y] = new_position(U,V,imax,jmax,delx,dely,delt,part_x,part_y)
    if( length(part_x) ~= length(part_y) )
        fprintf('Particle vectors differ in length!');
        return;
    end
    % Loop over all particles, compute individual u and v with 
    % interpolation() and propagate the particle
    for i = 1:length(part_x)
        [u,v] = interpolation(U,V,imax,jmax,delx,dely,part_x(i),part_y(i));
        part_x(i) = part_x(i) + u*delt;
        part_y(i) = part_y(i) + v*delt;
    end
end

