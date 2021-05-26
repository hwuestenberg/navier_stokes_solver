%######################################################
% Inject particles on the left-hand side of the flow. Creating new vectors
% Last: Created
%
% Created on 16/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [part_x,part_y] = set_particles(N, ug, og, delx)
    part_x = zeros(N,1);    % Pre-allocate memory
    part_y = zeros(N,1);
    
    for i = 1:N
        part_x(i) = delx;           % Avoid problems at boundary
        part_y(i) = (og-ug)/N*i+ug; % delta/N*i + offset, i = 1...N
    end
end

