Visualisation - Assignment 6


1)Particletrace.m

	In particletace.m first of all new_position function has been used , which calls interpolation function inside it for calculating
the velocities in the cell with the use of four corner values of a cell. In interpolation the removed particles have been ignored 
with setting velocities to 0.

	After interpolation , new position of particles are calculated with the use of Euler-equatiom.

	After new_position the particles leaving the grid are removed and according to the given problem , the particles intersecting with
obstacles are also removed.

	All the new postion of particles is stored in particle_trace matrix as columns.

2)Streaklines.m

	In streaklines calculation after each 4 time steps new particles are injected. So after each 4 time steps set_particles fumction 
is called for injection and each previou particles are moved one column ahead in particle_streak matrix.

	Meanwhile already injected particles are updated at every time step.As new particles are added these new particles will also b
be updated from that current time.

 
3)Particlevis.m

	This shows the visualisation of the pathlines and streaklines.

	Each path line is shown as the row of particle_x_trace  and particle_y_trace as x and y axis values.

	Same concept is used for streaklines with particle_x_streak and particle_y_streak matrices.





Multigrid Solver - Assignment 7


The multigrid solvers are called within the incompvis.m script. One can find them within the if-statements,

that evaluate the 'solver' parameter.


	- The mulitgrid solvers are called using a caller function multigrid.m. This function evaluates the parameter 'solver'.

	Depending on the chosen solver, the caller function calls a V-Cycle, W-Cycle or Full-Multigrid solver function.

	These functions call themselves recursively to run the given algorithms.


	- The smoothing iterations ('iter1' and 'iter2') are set at the top of the caller function multigrid.m.


Comparisons between different iterative schemes are made at the end of the incompvis.m script. 

	- The function 'evaluate_solver.m' plots the averaged times, iterations and residuals in a bar diagram. 
	
	Additionally, the convergence behaviour is plotted i.e. 2-norm of the residual over the iterations.

