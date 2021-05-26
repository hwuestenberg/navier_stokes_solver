#ifndef FLUID_HPP_
#define FLUID_HPP_

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::last;
using Eigen::all;
using Eigen::HouseholderQR;

class Fluid
{
public:
    int imax, jmax;
    double dx, dy, dt, alpha, nu, GX, GY, T_end;
    MatrixXd X, Y, U, V, P, F, G, K2D, RHS;
    Fluid(int problemID);
    ~Fluid();
    void setup();
    void run();
    void write_results();
};

// Setup
MatrixXd toeplitz(int imax, double factor);
MatrixXd kronecker_product(MatrixXd& A, MatrixXd& B);
double adjust_time_step_for_stability(double dx, double dy, double dt, double nu);

// Run
bool check_cfl_condition(MatrixXd &U, MatrixXd &V, double dx, double dy, double dt);
void set_boundary_conditions(MatrixXd &U, MatrixXd &V, MatrixXd &P, MatrixXd &F, MatrixXd &G);

void get_rhs(MatrixXd &U, MatrixXd &V, MatrixXd &F, MatrixXd &G, MatrixXd &RHS, int imax, int jmax, double dx, double dy, double dt, double GX, double GY, double alpha, double nu);
void assemble_rhs(MatrixXd &F, MatrixXd &G, MatrixXd &RHS, int imax, int jmax, double dx, double dy, double dt);
void assemble_fg(MatrixXd &U, MatrixXd &V, MatrixXd &F, MatrixXd &G, int imax, int jmax, double dx, double dy, double dt, double GX, double GY, double alpha, double nu);

void linear_solver(MatrixXd &K2D, MatrixXd &RHS, MatrixXd &P, int imax, int jmax);
void update_UV(MatrixXd &U, MatrixXd &V, MatrixXd &P, MatrixXd &F, MatrixXd &G, int imax, int jmax, double dx, double dy, double dt);

#endif