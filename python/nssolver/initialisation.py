import numpy as np
from scipy.linalg import toeplitz


def init_grid(imax, jmax):
    U = np.zeros((imax + 2, jmax + 2))
    V = np.zeros((imax + 2, jmax + 2))
    P = np.zeros((imax + 2, jmax + 2))
    F = np.zeros((imax + 2, jmax + 2))
    G = np.zeros((imax + 2, jmax + 2))

    return U, V, P, F, G


def init_linear(imax, jmax, dx, dy):
    columny = np.array([2, -1, *[0 for i in range(jmax-2)]])
    K1Dy = 1 / dy**2 * toeplitz(columny)    # Create Toeplitz matrix for one dimension
    Iy = np.identity(jmax)                # Identity matrix for block tridiagonal matrix

    columnx = np.array([2, -1, *[0 for i in range(imax-2)]])
    K1Dx = 1 / dx**2 * toeplitz(columnx)
    Ix = np.identity(imax)

    K2D = np.kron(K1Dy, Ix) + np.kron(Iy, K1Dx)  # Create block tridiagonal matrix for two dimensions

    return K2D


def get_problem_parameter(problemID):
    # 1 = Driven Cavity
    # 2 = Van Karman vortex street
    # 3 = Flow above a Stair
    # 4 = Tunnel with Constricitons
    if problemID == 1:
        problem_str = "Driven Cavity"

        # Discretisation
        imax = 16
        jmax = 16
        dx = 0.2
        dy = 0.2
        dt = 0.01
        T_end = 1

        # Solver
        alpha = 0.12

        # Fluid properties, initial velocity and pressure field
        nu = 0.4
        GX = 0
        GY = 0

        # Boundary conditions
        # 2 = No-slip
        wl = 2     # Left wall
        wr = 2     # Right wall
        wt = 2     # Top wall
        wb = 2     # Bottom wall

    elif problemID == 2:
        problem_str = "Von Karman Vortex Street"

        # Discretisation
        imax = 64
        jmax = 16
        dx = 0.2
        dy = 0.2
        dt = 0.02
        T_end= 20

        # Solver
        alpha = 0.9

        # Fluid properties, initial velocity and pressure field
        nu = 1.3414e-5
        GX = 0
        GY = 0

        # Boundary conditions
        # 1 = Free slip
        # 3 = Inflow/Outflow (Inflow at left wall)
        wl = 3     # Left wall
        wr = 3     # Right wall
        wt = 1     # Top wall
        wb = 1     # Bottom wall

    elif problemID == 3:
        problem_str = "Flow above a Stair"

        # Discretisation
        imax = 64
        jmax = 16
        dx = 0.2
        dy = 0.2
        dt = 0.02
        T_end= 20

        # Solver
        alpha = 0.12

        # Fluid properties, initial velocity and pressure field
        nu = 0.1
        GX = 0
        GY = 0

        # Boundary conditions
        # 2 = No-slip
        # 3 = Inflow/Outflow (Inflow at left wall)
        wl = 3     # Left wall
        wr = 3     # Right wall
        wt = 2     # Top wall
        wb = 2     # Bottom wall

    elif problemID == 4:
        problem_str = "Tunnel with Constrictions"

        # Discretisation
        imax = 64
        jmax = 16
        dx = 0.2
        dy = 0.2
        dt = 0.02
        T_end= 0.5

        # Solver
        alpha = 0.9

        # Fluid properties, initial velocity and pressure field
        nu = 1.3414e-5
        GX = 0
        GY = 0

        # Boundary Conditions
        # 1 = Free slip
        # 3 = Inflow/Outflow (Inflow at left wall)
        wl = 3     # Left wall
        wr = 3     # Right wall
        wt = 1     # Top wall
        wb = 1     # Bottom wall

    else:
        assert False, "Undefined problem ID!"

    return imax, jmax, dx, dy, dt, T_end, alpha, nu, GX, GY
