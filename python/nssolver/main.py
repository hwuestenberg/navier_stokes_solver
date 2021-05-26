#######################################################
# This script executes a Navier-Stokes-based simulation for incompressible,
# viscid, steady and two-dimensional flow.
#
# Created on 27/05/2018 with Matlab R2017b
# Ported to Python 3.8 on 03/05/2021
# by Henrik Wuestenberg
#######################################################
import numpy as np
import matplotlib.pyplot as plt
from time import time

from initialisation import *
from solver import *
from visual import *


if __name__ == '__main__':
    imax, jmax, dx, dy, dt, T_end, alpha, nu, GX, GY = get_problem_parameter(1)

    # Adjust time step for stable scheme
    check_stability(dx, dy, dt, nu)

    U, V, P, F, G = init_grid(imax, jmax)
    K2D = init_linear(imax, jmax, dx, dy)
    # fig, ax = init_plot(U, V, imax, jmax)

    T = 0               # Physical time
    t_start = time()    # Computational time
    while T < T_end:
        # Set boundary conditions
        U, V, P, F, G = boundary_conditions(U, V, P, F, G)

        # Solve pressure equation (w/ Gaussian Elimination)
        RHS = get_RHS(U, V, F, G, imax, jmax, dx, dy, dt, GX, GY, alpha, nu)
        P = linear_solver(K2D, RHS, P, imax, jmax)
        U, V = update_UV(U, V, P, F, G, imax, jmax, dx, dy, dt)

        # Check stability and, if necessary, stop simulation
        assert cfl_condition(U, V, dx, dy, dt), "CFL condition not satisfied"

        # update_plot(U, V, imax, jmax, ax)  # Animated simulation

        # Time stepping
        T = T + dt
        # print("Physical time: {t_phy}".format(t_phy=T))

    t_total = time() - t_start
    print("Computational Time: {t} seconds".format(t=np.round(t_total)))
    # update_plot(U, V, imax, jmax, T, ax)  # Plot with streamlines
    fig, ax = init_plot(U, V, imax, jmax)
    plt.pause(5)

    # Compare to cpp results
    Xpp, Ypp, Upp, Vpp, Ppp = load_cpp_results(imax, jmax)
    figpp, axpp = init_plot(Upp, Vpp, imax, jmax)
    plt.pause(5)

