#include "fluid.hpp"


Fluid::Fluid(int problemID)
{
    // Problem setup
    switch (problemID)
    {
    case 1:
        // Physical domain
        imax = 16;
        jmax = 16;
        dx = 0.2;
        dy = 0.2;
        dt = 0.01;
        T_end = 1;

        // Solver
        alpha = 0.12;

        // Fluid properties, initial velocity and pressure field
        nu = 0.4;
        GX = 0;
        GY = 0;
        break;
    
    default:
        std::cout << "Invalid problemID: " << problemID << std::endl;
        assert(0);
        break;
    }
}

Fluid::~Fluid()
{
}

void Fluid::setup() 
{
    std::cout << "Setting up matrices" << std::endl;

    // Check stability of numerical scheme
    assert(adjust_time_step_for_stability(dx, dy, dt, nu));

    // Initialise grid
    VectorXd x = VectorXd::LinSpaced(imax + 2, 0.0, dx * (imax + 2));
    VectorXd y = VectorXd::LinSpaced(jmax + 2, 0.0, dy * (jmax + 2));
    // std::cout << "x" << x.rows() << x.cols() << std::endl << x << std::endl << std::endl;
    // std::cout << "y" << y.rows() << y.cols() << std::endl << y << std::endl << std::endl;
    X = x.replicate(jmax + 2, 1).reshaped(imax + 2, jmax + 2);
    Y = y.replicate(imax + 2, 1).reshaped(imax + 2, jmax + 2).transpose();
    // std::cout << "X" << X.rows() << X.cols() << std::endl << X << std::endl << std::endl;
    // std::cout << "Y" << Y.rows() << Y.cols() << std::endl << Y << std::endl << std::endl;

    U = MatrixXd::Zero(imax + 2, jmax + 2);
    V = MatrixXd::Zero(imax + 2, jmax + 2);
    P = MatrixXd::Zero(imax + 2, jmax + 2);
    F = MatrixXd::Zero(imax + 2, jmax + 2);
    G = MatrixXd::Zero(imax + 2, jmax + 2);

    // Initialise linear system (Block matrix)
    MatrixXd K1Dx = toeplitz(imax, 1 / pow(dx, 2.0));
    MatrixXd K1Dy = toeplitz(jmax, 1 / pow(dy, 2.0));
    // std::cout << K1Dx << std::endl << std::endl;
    // std::cout << K1Dy << std::endl;
    
    MatrixXd Ix = MatrixXd::Identity(imax, imax);
    MatrixXd Iy = MatrixXd::Identity(jmax, jmax);
    // std::cout << Ix << std::endl << std::endl;
    // std::cout << Iy << std::endl;

    K2D = kronecker_product(K1Dy, Ix) + kronecker_product(Iy, K1Dx);
    // std::cout << "K2D shape NxM: " << K2D.rows() << "x" << K2D.cols() << std::endl;
    // std::cout << K2D << std::endl;

    std::cout << "Setup completed" << std::endl;
}

void Fluid::run() 
{
    std::cout << "Code is running" << std::endl;

    double T = 0;                               // Physical time
    auto tStart = high_resolution_clock::now(); // Computational time

    while(T < T_end) 
    {
        std::cout << "Physical time: " << std::to_string(T) << std::endl;
        
        // Set boundary conditions
        set_boundary_conditions(U, V, P, F, G);

        // // Solve pressure equation (w/ Gaussian Elimination)
        get_rhs(U, V, F, G, RHS, imax, jmax, dx, dy, dt, GX, GY, alpha, nu);
        linear_solver(K2D, RHS, P, imax, jmax);
        update_UV(U, V, P, F, G, imax, jmax, dx, dy, dt);

        // Check stability and, if necessary, stop simulation
        assert(check_cfl_condition(U, V, dx, dy, dt));

        // Time stepping
        T = T + dt;
    }
    auto tEnd = high_resolution_clock::now();
    std::cout << "Simulation finished at physical time: " << std::to_string(T) << std::endl;
    
    duration<double> tComp = duration_cast<duration<double>>(tEnd - tStart);
    std::cout << "Computational time: " << tComp.count() << " s" << std::endl;
}

void Fluid::write_results() 
{
    std::cout << "Writing results to file results_cpp.csv" << std::endl;

    std::ofstream file("./results_cpp.csv", std::ofstream::out);
    if(file.is_open())
    {
        // Header
        file << "X,Y,U,V,P" << std::endl;
        
        // Data
        for( int i = 0; i < imax + 2; i++)
        {
            for( int j = 0; j < jmax + 2; j++)
            {
                file << X(i, j) << "," << Y(i, j) << "," << U(i, j) << "," << V(i, j) << "," << P(i, j) << std::endl;
            }
        }
    }
    file.close();

    std::cout << "Results saved" << std::endl;
}

MatrixXd toeplitz(int size, double factor) 
{
    MatrixXd A = MatrixXd::Zero(size, size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if ( j == i )
            {
                A(i, j) = 2 * factor;
            }
            else if ( j == i + 1 || j == i - 1)
            {
                A(i, j) = -1 * factor;
            }
        }
        
    }
    return A;
}

MatrixXd kronecker_product(MatrixXd &A, MatrixXd &B)
{
    int Arows = A.rows();
    int Acols = A.cols();
    int Brows = B.rows();
    int Bcols = B.cols();

    // std::cout << "Ai " << Arows << "  Aj " << Acols << std::endl;
    // std::cout << "Bk " << Brows << "  Bl " << Bcols << std::endl;

    MatrixXd C = MatrixXd::Zero(Arows * Brows, Acols * Bcols);

    for (int i = 0; i < Arows; i++)
    {
        for (int k = 0; k < Brows; k++)
        {
            for (int j = 0; j < Acols; j++)
            {
                for (int l = 0; l < Bcols; l++)
                {
                    // std::cout << i << "" << j << "" << k << "" << l << "  " << k + i * Brows << "" << l + j * Bcols << std::endl;
                    C(k + i * Brows, l + j * Bcols) = A(i, j) * B(k, l);
                }
                
            }
        }
        
    }
    return C;
    
}

double adjust_time_step_for_stability(double dx, double dy, double dt, double nu)
{
    bool condition_value;
    // Evaluate stability criteria once
    condition_value = 2 * nu * dt < pow(dy, 2.0) / (1 + pow(dy / dx, 2.0));
    std::cout << "Initial condition value is " << condition_value << std::endl;

    // Decrease time step size if stability condition not satisfied
    while( !condition_value )
    {
        dt = dt - dt / 2;
        condition_value = 2 * nu * dt < pow(dy, 2.0) / (1 + pow(dy / dx, 2.0));
        std::cout << "time step adjusted to " << dt << std::endl;
    }
    std::cout << "After adjustment; condition value is " << condition_value << std::endl;

    // Output resulting time step size or error
    assert(condition_value);

    return dt;
}

bool check_cfl_condition(MatrixXd &U, MatrixXd &V, double dx, double dy, double dt)
{
    bool cfl_x, cfl_y;
    // std::cout << "cwiseAbs: " << std::endl << U.cwiseAbs() << std::endl;
    // std::cout << "cwiseAbs.maxCoeff: " << std::endl << U.cwiseAbs().maxCoeff() << std::endl;
    // std::cout << "cwiseAbs.maxCoeff * dt: " << std::endl << U.cwiseAbs().maxCoeff() * dt << std::endl;
    cfl_x = U.cwiseAbs().maxCoeff() * dt < dx;
    cfl_y = V.cwiseAbs().maxCoeff() * dt < dy;
    return true;
}

void set_boundary_conditions(MatrixXd &U, MatrixXd &V, MatrixXd &P, MatrixXd &F, MatrixXd &G)
{
    // Matrices for convenience
    double lidVel = 2;
    MatrixXd constVel = MatrixXd::Constant(U.rows() - 2, 1, lidVel);  // Upper wall velocity
    MatrixXd zCols = MatrixXd::Zero(1, U.cols());    // Zero velocity
    // std::cout << constVel.size() << std::endl;

    // No-slip on U
    U(all, 0) = -U(all, 1); // Lower boundary
    U(seq(1, last - 1), last) = constVel - U(seq(1, last - 1), last - 1);  // Upper boundary
    U(0, all) = zCols;  // Left boundary
    U(last - 1, all) = zCols;  // Right boundary
    // std::cout << U << std::endl;

    // No-slip on V
    V(all, 0) = zCols.transpose();  // Lower boundary
    V(all, last - 1) = zCols.transpose();  // Upper boundary
    V(0, all) = -V(1, all);  // Left boundary
    V(last, all) = -V(last - 1, all);  // Right boundary
    // std::cout << V << std::endl;

    // No-slip on P
    P(all, 0) = P(all, 1);  // Lower boundary
    P(all, last) = P(all, last - 1);  // Upper boundary
    P(0, all) = P(1, all);  // Left boundary
    P(last, all) = P(last - 1, all);  // Right boundary
    // std::cout << P << std::endl;

    // BC for virtual term F
    F(0, all) = U(0, all);
    F(last - 1, all) = U(last - 1, all);
    // std::cout << F << std::endl;

    // BC for virtual term G
    G(all, 0) = V(all, 0);
    G(all, last -1) = V(all, last -1);
    // std::cout << G << std::endl;
}

void get_rhs(MatrixXd &U, MatrixXd &V, MatrixXd &F, MatrixXd &G, MatrixXd &RHS, int imax, int jmax, double dx, double dy, double dt, double GX, double GY, double alpha, double nu)
{
    assemble_fg(U, V, F, G, imax, jmax, dx, dy, dt, GX, GY, alpha, nu);
    // std::cout << F << std::endl << std::endl;
    // std::cout << G << std::endl << std::endl;

    assemble_rhs(F, G, RHS, imax, jmax, dx, dy, dt);
    // std::cout << RHS << std::endl;
}

void assemble_rhs(MatrixXd &F, MatrixXd &G, MatrixXd &RHS, int imax, int jmax, double dx, double dy, double dt)
{
    double f, g;    // For convenience
    RHS = MatrixXd::Zero(imax + 2, jmax + 2);

    for( int i = 1; i < imax + 1; i++)
    {
        for( int j = 1; j < jmax + 1; j++)
        {
            f = F(i - 1, j) - F(i, j);
            g = G(i, j - 1) - G(i, j);
            RHS(i, j) = (f / dx + g / dy) / dt;
        }
    }

}

void assemble_fg(MatrixXd &U, MatrixXd &V, MatrixXd &F, MatrixXd &G, int imax, int jmax, double dx, double dy, double dt, double GX, double GY, double alpha, double nu)
{
    double d2u, duu_dx, duv_dy, fbracket;
    double d2v, dvv_dy, duv_dx, gbracket;
    
    for( int i = 1; i < imax; i++)
    {
        for( int j = 1; j < jmax + 1; j++)
        {
            d2u = (U(i + 1, j) - 2 * U(i, j) + U(i - 1, j)) / dx / dx + (U(i, j + 1) - 2 * U(i, j) + U(i, j - 1)) / dy / dy;
            duu_dx = (pow((U(i, j) + U(i + 1, j)) / 2, 2.0) - pow((U(i - 1, j) + U(i, j)) / 2, 2.0)) / dx + 
            alpha / dx * (abs(U(i, j) + U(i + 1, j)) / 2 * (U(i, j) - U(i + 1, j)) / 2 - abs(U(i - 1, j) + U(i, j)) / 2 * (U(i - 1, j) - U(i, j)) / 2);
            duv_dy = ((V(i, j) + V(i + 1, j)) / 2 * (U(i, j) + U(i, j + 1)) / 2 - (V(i, j - 1) + V(i + 1, j - 1)) / 2 * (U(i, j - 1) + U(i, j)) / 2) / dy + 
            alpha / dy * (abs(V(i, j) + V(i + 1, j)) / 2 * (U(i, j) - U(i, j + 1)) / 2 - abs(V(i, j - 1) + V(i + 1, j - 1)) / 2 * (U(i, j - 1) - U(i, j)) / 2);
            fbracket = nu * d2u - duu_dx - duv_dy + GX;
            F(i, j) = U(i, j) + dt * fbracket;
        }
    }

    for( int i = 1; i < imax + 1; i++)
    {
        for( int j = 1; j < jmax; j++)
        {
            d2v = (V(i + 1, j) - 2 * V(i, j) + V(i - 1, j)) / dx / dx + (V(i, j + 1) - 2 * V(i, j) + V(i, j - 1)) / dy / dy;
            dvv_dy = (pow((V(i, j) + V(i, j + 1)) / 2, 2.0) - pow((V(i, j - 1) + V(i, j)) / 2, 2.0)) / dy + 
            alpha / dy * (abs(V(i, j) + V(i, j + 1)) / 2 * (V(i, j) - V(i, j + 1)) / 2 - abs(V(i, j - 1) + V(i, j)) / 2 * (V(i, j - 1) - V(i, j)) / 2);
            duv_dx = ((U(i, j) + U(i, j + 1)) / 2 * (V(i, j) + V(i + 1, j)) / 2 - (U(i - 1, j) + U(i - 1, j + 1)) / 2 * (V(i - 1, j) + V(i, j)) / 2) / dx + 
            alpha / dx * (abs(U(i, j) + U(i, j + 1)) / 2 * (V(i, j) - V(i + 1, j)) / 2 - abs(U(i - 1, j) + U(i - 1, j + 1)) / 2 * (V(i - 1, j) - V(i, j)) / 2);
            gbracket = nu * d2v - dvv_dy - duv_dx + GY;
            G(i, j) = V(i, j) + dt * gbracket;
        }
    }
}

void linear_solver(MatrixXd &K2D, MatrixXd &RHS, MatrixXd &P, int imax, int jmax)
{
    double residual;
    MatrixXd rhs, p;

    rhs = RHS(seq(1, imax), seq(1, jmax)).reshaped(imax * jmax, 1);
    // std::cout << "Size rhs:\t" << rhs.size() << std::endl;
    // std::cout << "Size K2D:\t" << K2D.size() << std::endl;

    p = K2D.householderQr().solve(rhs);     // QR + forward substitution (Direct solver)
    
    residual = (K2D * p - rhs).norm() / rhs.norm();     // norm == L2-norm
    assert(residual <= 1e-5);

    // std::cout << P << std::endl << std::endl;
    P(seq(1, imax), seq(1, jmax)) = p.reshaped(imax, jmax);
    // std::cout << P << std::endl << std::endl;
}

void update_UV(MatrixXd &U, MatrixXd &V, MatrixXd &P, MatrixXd &F, MatrixXd &G, int imax, int jmax, double dx, double dy, double dt)
{
    // Update U
    for( int i = 1; i < imax; i++)
    {
        for( int j = 1; j < jmax + 1; j++)
        {
            U(i, j) = F(i, j) - dt / dx * (P(i + 1, j) - P(i, j));
        }
    }

    // Update V
    for( int i = 1; i < imax + 1; i++)
    {
        for( int j = 1; j < jmax; j++)
        {
            V(i, j) = G(i, j) - dt / dy * (P(i, j + 1) - P(i, j));
        }
    }
}


/*
def init_plot(U, V, imax, jmax):
    plt.ion()
    figure = plt.figure()
    axis = figure.subplots()
    axis.contourf(np.arange(imax + 2), np.arange(jmax + 2), np.sqrt(U.T**2 + V.T**2))
    axis.set_aspect("equal")
    return figure, axis


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
*/