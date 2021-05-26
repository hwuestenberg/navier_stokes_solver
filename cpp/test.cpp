// Created on 08/05/2021 by hw
#include <iostream>
#include <Eigen/Dense>
#include "fluid.hpp"

using Eigen::MatrixXd;
using Eigen::seq;
using Eigen::last;
using Eigen::all;

// int arr[10];
// int* ptr = new int[10];

// class example
// {
// public:
//     int a, b, c;
//     example(int value) {
//         a = value;
//         b = 2;
//     }
//     ~example() {};
//     void do_cout()
//     {
//         std::cout << "[a, b] = [" << std::to_string(a) << ", " << std::to_string(b) << "]" << std::endl;
//     }
//     void change_b(int value)
//     {
//         b = value;
//     }
// };

// void use_reference(example& ex)
// {
//     ex.change_b(-4);
// }

// void use_copy(example ex)
// {
//     ex.change_b(0);
// }

int main() {
    std::cout << "Hello World!\n";
    

    // // Array stuff
    // for (int i = 0; i < 10; i++) {
    //     arr[i] = i;
    //     *(ptr + i) = i;// Equivalent to ptr[i]
    // }

    // for (int j = 0; j < 10; j++) {
    //     std::cout << arr[j];
    // }
    // std::cout << std::endl;
    
    // for (int k = 0; k < 10; k++) {
    //     std::cout << ptr[k];
    // }
    // std::cout << std::endl;
    

    // // Pointer stuff
    // int j = 2;
    // std::cout << " j      " << j << std::endl;  // Integer
    // // std::cout << "*j      " << *j << std::endl;  Throws error > not a pointer
    // std::cout << "&j      " << &j << std::endl;  // Address of integer

    // // int* ptri = new int(1);
    // int* ptri = &j;
    // std::cout << " ptri   " << ptri << std::endl;  // Address of integer
    // std::cout << "&ptri   " << &ptri << std::endl; // Address of pointer
    // std::cout << "*ptri   " << *ptri << std::endl;  // Integer
    // std::cout << "&*ptri   " << &(*ptri) << std::endl;  // Address of integer

    // example* e = new example(3);
    // e->do_cout();

    // use_reference(*e);
    // (*e).do_cout();  // (*e).f() == e->f()
    
    // use_copy(*e);
    // e->do_cout();


    // Equivalent loops
    // Python   for i in range(1, 11)
    // C++      for(i = 1; i < 11; i++)


    // Non-equivalent slicing
    // Python   1:N
    // C++      seq(1, N-1)


// TODO CORRECT indices for F, G in Python and Matlab version


    // Eigen/Dense stuff
    // int N = 4, M = 4;
    // MatrixXd m(N, M), ms(N-2, M-2);
    // for (int i = 0; i < N; i++) {
    //     for ( int j = 0; j < M; j++) {
    //         m(i, j) = i * j;
    //     }
    // }
    // std::cout << m << std::endl << std::endl;
    // std::cout << m.reshaped(M*N, 1) << std::endl << std::endl;
    // std::cout << m(seq(1, last), 1) << std::endl;
    // std::cout << m(seq(1, N - 2), seq(1, M - 2)) << std::endl << std::endl;

    // ms = m(seq(1, N - 2), seq(1, M - 2)).reshaped(2, 2);
    // // do stuff with sub-matrix
    // ms << 7, 7,
    //       7, 7;
    // std::cout << ms << std::endl << std::endl;

    // m(seq(1, N - 2), seq(1, M - 2)) = ms;
    // std::cout << m << std::endl << std::endl;

    // MatrixXd b(2, 1), A(2, 2), x;
    // b << 1, 2;
    // A << 1, 2,
    //      3, 5;
    // x = A.householderQr().solve(b);
    // std::cout << x << std::endl;


    // // Class Fluid (compile with "g++ test.cpp fluid.cpp -o test")
    Fluid fluid(1);
    fluid.setup();
    fluid.run();
    fluid.write_results();
    // // assert(0);

    return 0;
}