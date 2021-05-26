%######################################################
% Check stability condition 26.1 and, eventually, decrease time step size.
%
% Created on 27/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function stab1(delt,delx,dely,nu)
    % Evaluate stability criteria once
    stab1 = 2*nu*delt < dely^2/(1+(dely/delx)^2);
    
    % Iteratively decrease time step size if stability condition does not 
    % hold
    while stab1 ~= 1
        delt = delt - delt/2;
        stab1 = 2*nu*delt < dely^2/(1+(dely/delx)^2);
    end
    
    % Output resulting time step size or error
    if stab1 ~= 1
        disp('Stability conidition 26.1 is not fulfilled for any time step size');
    else
        stab1_str = sprintf('The time step is determined to %5.3f for stability',delt);
        disp(stab1_str);
    end
end