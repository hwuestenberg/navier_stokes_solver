# navier_stokes_solver
Finite differences solver for 2D Navier-Stokes equations

Model assumptions:
- incompressible,
- unsteady,
- two-dimensional.



Tested several solver for linear systems of the form Ax=b, including:
- Successive-over-Relaxation,
- Conjugate Gradients and
- Mutligrid (V- and W-cycle)

Applied to canonical flows, including:
- Backwards-facing step,
- Sudden Converging-Diverging Channel,
- Flow around a Square Cylinder and
- Lid-Driven Cavity
