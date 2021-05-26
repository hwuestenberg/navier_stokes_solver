import numpy as np
from scipy.linalg import solve, norm


def linear_solver(K2D, RHS, P, imax, jmax):
    rhs = np.reshape(RHS[1:imax + 1, 1:jmax + 1], (imax * jmax, 1))

    p = solve(K2D, rhs)
    res = norm(np.dot(K2D, p) - rhs, 2)
    assert res <= 1e-5, "Residual of direct solver is larger than 1e-5!"

    P[1:imax + 1, 1:jmax + 1] = np.reshape(p, (imax, jmax))

    return P


def get_RHS(U, V, F, G, imax, jmax, dx, dy, dt, GX, GY, alpha, nu):
    F, G = calc_FG(U, V, F, G, imax, jmax, dx, dy, dt, GX, GY, alpha, nu)
    RHS = calc_RHS(F, G, imax, jmax, dx, dy, dt)
    return RHS


def calc_RHS(F, G, imax, jmax, dx, dy, dt):
    RHS = np.zeros((imax + 2, jmax + 2))

    for i in range(1, imax + 1):
        for j in range(1, jmax + 1):
            f = F[i - 1, j] - F[i, j]
            g = G[i, j - 1] - G[i, j]
            RHS[i, j] = (f / dx + g / dy) / dt
    return RHS


def calc_FG(U, V, F, G, imax, jmax, dx, dy, dt, GX, GY, alpha, nu):
    for i in range(1, imax):
        for j in range(1, jmax+1):
            d2u = (U[i + 1, j] - 2 * U[i, j] + U[i - 1, j]) / dx / dx + (U[i, j + 1] - 2 * U[i, j] + U[i, j - 1]) / dy / dy
            duu_dx = (((U[i, j] + U[i + 1, j]) / 2) ** 2 - ((U[i - 1, j] + U[i, j]) / 2) ** 2) / dx + \
                     alpha / dx * ((abs(U[i, j] + U[i + 1, j])) / 2 * (U[i, j] - U[i + 1, j]) / 2 - (
                             abs(U[i - 1, j] + U[i, j])) / 2 * (U[i - 1, j] - U[i, j]) / 2)
            duv_dy = ((V[i, j] + V[i + 1, j]) / 2 * (U[i, j] + U[i, j + 1]) / 2 - (V[i, j - 1] + V[i + 1, j - 1]) / 2 * (U[i, j - 1] + U[i, j]) / 2) / dy + \
                     alpha / dy * ((abs(V[i, j] + V[i + 1, j])) / 2 * (U[i, j] - U[i, j + 1]) / 2 - (
                             abs(V[i, j - 1] + V[i + 1, j - 1])) / 2 * (U[i, j - 1] - U[i, j]) / 2)
            Fbracket = nu * d2u - duu_dx - duv_dy + GX
            F[i, j] = U[i, j] + dt * Fbracket

    for i in range(1, imax+1):
        for j in range(1, jmax):
            d2v = (V[i + 1, j] - 2 * V[i, j] + V[i - 1, j]) / dx / dx + (V[i, j + 1] - 2 * V[i, j] + V[i, j - 1]) / dy / dy
            dvv_dy = (((V[i, j] + V[i, j + 1]) / 2) ** 2 - ((V[i, j - 1] + V[i, j]) / 2) ** 2) / dy + \
                     alpha / dy * ((abs(V[i, j] + V[i, j + 1])) / 2 * (V[i, j] - V[i, j + 1]) / 2 - (
                    abs(V[i, j - 1] + V[i, j])) / 2 * (V[i, j - 1] - V[i, j]) / 2)
            duv_dx = ((U[i, j] + U[i, j + 1]) / 2 * (V[i, j] + V[i + 1, j]) / 2 - (U[i - 1, j] + U[i - 1, j + 1]) / 2 * (V[i - 1, j] + V[i, j]) / 2) / dx + \
                     alpha / dx * (
                             (abs(U[i, j] + U[i, j + 1])) / 2 * (V[i, j] - V[i + 1, j]) / 2 -
                             (abs(U[i - 1, j] + U[i - 1, j + 1])) / 2 * (V[i - 1, j] - V[i, j]) / 2)
            Gbracket = nu * d2v - dvv_dy - duv_dx + GY
            G[i, j] = V[i, j] + dt * Gbracket

    return F, G


def update_UV(U, V, P, F, G, imax, jmax, dx, dy, dt):
    for i in range(1, imax):
        for j in range(1, jmax + 1):
            U[i, j] = F[i, j] - dt / dx * (P[i + 1, j] - P[i, j])

    for i in range(1, imax + 1):
        for j in range(1, jmax):
            V[i, j] = G[i, j] - dt / dy * (P[i, j + 1] - P[i, j])

    return U, V


def boundary_conditions(U, V, P, F, G):
    # No-slip on U
    U[:, 0] = -U[:, 1]  # Lower boundary
    U[1:-1, -1] = 2 - U[1:-1, -2]  # Upper boundary
    U[0, :] = 0  # Left boundary
    U[-2, :] = 0  # Right boundary

    # No-slip on V
    V[:, 0] = 0  # Lower boundary
    V[:, -2] = 0  # Upper boundary
    V[0, :] = -V[1, :]  # Left boundary
    V[-1, :] = -V[-2, :]  # Right boundary

    # No-slip on P
    P[:, 0] = P[:, 1]  # Lower boundary
    P[:, -1] = P[:, -2]  # Upper boundary
    P[0, :] = P[1, :]  # Left boundary
    P[-1, :] = P[-2, :]  # Right boundary

    # BC for virtual term F
    F[0, :] = U[0, :]
    F[-2, :] = U[-2, :]

    # BC for virtual term G
    G[:, 0] = V[:, 0]
    G[:, -2] = V[:, -2]

    return U, V, P, F, G


def check_stability(dx, dy, dt, nu):
    # Evaluate stability criteria once
    scondition = 2 * nu * dt < dy ** 2 / (1 + (dy / dx) ** 2)

    # Decrease time step size if stability condition not satisfied
    while scondition is False:
        dt = dt - dt / 2
        scondition = 2 * nu * dt < dy ** 2 / (1 + (dy / dx) ** 2)

    # Output resulting time step size or error
    assert scondition, "Numerical scheme is unconditionally unstable!"

    return dt


def cfl_condition(U, V, dx, dy, dt):
    cfl_x = np.max(np.abs(U))*dt < dx
    cfl_y = np.max(np.abs(V))*dt < dy
    return cfl_x and cfl_y
