%######################################################
% Restriction operator for the multigrid method. Computes a coarse grid by 
% averaging four adjacent cells to a single cell. Level parameter 'L'
% recognises the depth of current grid. L_0 = 0
% Works only for rectangular grids with imax > jmax
%
%
% Last: Created
%
%
% Created on 10/08/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function coarse = restriction(fine,imax,jmax,L)
    L = 2^L;
    new_size = length(fine)/4;
    coarse = zeros(new_size,1);
    
    y_offset = (imax-1)/L;  % Offset to j+1 values in rhs-vector
    x_offset = 0;           % Averaging offset, to skip elements of j = 2,4,6,...
    
    for i=1:new_size        % Construct average from four adjacent cells
        a = 2*i+x_offset;   % Coefficient for convenience
        
        % Calculate average using cells: (i,j),(i+1,j),(i,j+1),(i+1,j+1)
        if (jmax-1)/L ~= 1
            coarse(i) = ( fine(a-1)+fine(a)+fine(a-1+y_offset)+fine(a+y_offset) )/4;
        
        % Use just 2 cells, if smaller dimension has length=1
        else
            coarse(i) = ( fine(a-1)+fine(a) )/2;
            
        end
        if( a == y_offset+2*x_offset ) % Check end of row
            x_offset = x_offset + y_offset;
        end
    end
end

