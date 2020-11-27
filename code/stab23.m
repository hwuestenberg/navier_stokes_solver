%######################################################
% Check CFL condition for stability i.e. 26.2 and 26.3
%
% Created on 29/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function ret = stab23(U,V,delt,delx,dely,T)
    % Evaluate stability criteria once
    stab2 = max(max(abs(U)))*delt < delx;
    stab3 = max(max(abs(V)))*delt < dely;
        
    % Output if simulation becomes unstable
    if ~stab2
        fprintf('StabCond 26.2 not fulfilled at %4.3f seconds!\n',T);
        ret = 0;
    else
        ret = 1;
    end
    
    if ~stab3
        fprintf('StabCond 26.3 not fulfilled at %4.3f seconds!\n',T);
        ret = 0;
    else
        ret = 1;
    end
end