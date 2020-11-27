%######################################################
% Calculate velocities for the next time step
% Last: Changed iteration starting from 2 instead of 1
%
% Created on 26/05/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function [U,V] = calc_UV(U,V,P,F,G,imax,jmax,delt,delx,dely)
    % Compute U
	for i = 2:imax-1
		for j = 2:jmax
			U(i,j) = F(i,j)-(delt/delx)*( P(i+1,j)-P(i,j) );
		end
	end
	% Compute V
	for i = 2:imax
		for j = 2:jmax-1
			V(i,j) = G(i,j)-(delt/dely)*( P(i,j+1)-P(i,j) );
		end
	end
	
end