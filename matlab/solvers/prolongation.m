%######################################################
% Prolongation operator for the multigrid method. The value from the corase
% grid is simply transferred to each daughter cell. Level parameter 'L'
% recognises the depth of current grid
%
%
% Last: Created
%
%
% Created on 11/08/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function fine = prolongation(coarse,imax,L)
    L = 2^L;                % Level parameter is power of 2
    new_size = length(coarse)*4;
    fine = zeros(new_size,1);
    
    y_offset = (imax-1)/L;  % Offset to j+1 values in rhs-vector
    x_offset = 0;           % Averaging offset, to avoid overlapping
    
    for i=1:length(coarse)  % Construct average from four adjacent cells
        a = 2*i+x_offset;   % Coefficient for convenience
        
        % Transfer coarse value to daughter cells
        fine(a-1)   = coarse(i);
        fine(a)     = coarse(i);
        fine(a-1+y_offset)  = coarse(i);
        fine(a+y_offset)    = coarse(i);
        
        if( a == y_offset+2*x_offset ) % Check end of row
            x_offset = x_offset + y_offset;
        end
    end
end

