%######################################################
% Call the streakline subroutine ,injects particles after 
% each 4 time steps and position of all the particles is
% traced at each time step
% Last: Set removing conditions for problem 2 (16/07/2018)
%
% Created on 12/08/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################



function [ part_x_streak,part_y_streak ] = streaklines(U,V,imax,jmax,delx,dely,delt,part_x_streak,part_y_streak,N,ug,og,delt_n,T,T_end,problem )

%#----------------injecting new particles--------------#%
T=round(T,2);
%injects particles after each 4 time steps (delt_n*delt)

if((rem(T,delt_n*delt))<=10^(-10))
    for j=2:T_end/(delt_n*delt)
       k=T_end/(delt_n*delt)+2-j;
       part_x_streak(:,k)=part_x_streak(:,k-1);
       part_y_streak(:,k)=part_y_streak(:,k-1);
       
    end 
%part_x_streak(:,1)=0;
%part_y_streak(:,1)=0;    
[part_x_streak(:,1),part_y_streak(:,1)] = set_particles(N, ug, og, delx);
end

%#---------------updating already injected particles-----------#%

for i=1:floor(T/(delt_n*delt))+1
    
    [part_x_streak(:,i),part_y_streak(:,i)] = particletrace(U,V,imax,jmax,delx,dely,delt,part_x_streak(:,i),part_y_streak(:,i),problem);
    
end

