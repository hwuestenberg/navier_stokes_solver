%######################################################
% Plot Particle Tracer positions at every time step and Streaklines
% Last: Set removing conditions for problem 1&2 (16/07/2018)
%
% Created on 16/07/2018 with Matlab R2017b
% by Henrik Wüstenberg,Rushit Kansara, Sarvesh Vishishth
%######################################################
function evaluate_Solver(solv,epsi,omg,solver)
    
    cd ./solver/
    
    % Find current solver for title
    switch(solver)
        case 0
            str = 'Grid_Oriented';
        case 1
            str = 'SOR';
        case 2
            str = 'PCG';
        case 3
            str = 'Simple Multigrid';
        case 4
            str = 'V-Cycle Multigrid';
        case 5
            str = 'W-Cycle Multigrid';
        case 6
            str = 'Full Multigrid';
        case 7
            str = 'Direct';
    end
    
    % Create figure
    figure(12)
    title(sprintf('%s Solver, Averaged data',str));
    hold on
    
    % Calculate averags
    avg_time = mean(solv.time);
    avg_iter = mean(solv.iter);
    avg_res = mean(solv.res);
    averages = [avg_time avg_iter avg_res];
    
    % Plot averages
    names = categorical({'Time [s]','Iterations [-]','Residual [-]'});
    b = bar(names,averages);
    
    % Set colours
    b.FaceColor = 'flat';
    b.CData(1,:) = [0 0 1];
    b.CData(2,:) = [1 0 0];
    b.CData(3,:) = [0 1 0];
    
    set(gca, 'YScale', 'log')   % Change y-axis to logarithm
    set(gca, 'DefaultAxesFontSize',25)
    hold off
    
    print(sprintf('%s_epsi_%s_omg_%s_bar',str,num2str(log10(epsi)),num2str(100*omg)),'-dpng');
    
    % Plot convergence behaviour
    if ~isempty(find(solv.resvec,1))
        figure(13)
        semilogy(1:solv.iter(length(solv.iter)),solv.resvec(1:solv.iter(length(solv.iter))),'-or')
        
        % Set axis parameters
        title(sprintf('%s Solver',str));
        xlabel('Iteration')
        ylabel('2-Norm Residual')
        xlim([0 100])
        ylim([epsi*1e-1 1e+2])
        set(gca, 'DefaultAxesFontSize',25)
        
        print(sprintf('%s_epsi_%s_omg_%s_convergence',str,num2str(log10(epsi)),num2str(100*omg)),'-dpng');
    end
    
    
end