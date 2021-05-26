import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def init_plot(U, V, imax, jmax):
    plt.ion()
    figure = plt.figure()
    axis = figure.subplots()
    axis.contourf(np.arange(imax + 2), np.arange(jmax + 2), np.sqrt(U.T**2 + V.T**2))
    axis.set_aspect("equal")
    return figure, axis


def update_plot(U, V, imax, jmax, T, axis):
    axis.cla()
    X, Y = np.arange(imax + 2)[1:-1], np.arange(jmax + 2)[1:-1]
    axis.contourf(X, Y, np.sqrt(U[1:-1, 1:-1].T ** 2 + V[1:-1, 1:-1].T ** 2))
    plt.quiver(X, Y, U[1:-1, 1:-1].T, V[1:-1, 1:-1].T, color="w")
    axis.set_xlim([1, imax])
    axis.set_ylim([1, jmax])
    axis.text(imax + 1, jmax, "time {t:.3f}".format(t=T), bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    plt.draw()
    plt.pause(0.01)


def load_cpp_results(imax, jmax):
    df = pd.read_csv("../../cpp/results_cpp.csv")
    Xpp = df['X'].to_numpy().reshape(imax + 2, jmax + 2)
    Ypp = df['Y'].to_numpy().reshape(imax + 2, jmax + 2)
    Upp = df['U'].to_numpy().reshape(imax + 2, jmax + 2)
    Vpp = df['V'].to_numpy().reshape(imax + 2, jmax + 2)
    Ppp = df['P'].to_numpy().reshape(imax + 2, jmax + 2)

    return Xpp, Ypp, Upp, Vpp, Ppp
