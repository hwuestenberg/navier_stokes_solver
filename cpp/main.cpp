// Created on 08/05/2021 by hw
#include <iostream>
#include "fluid.hpp"

int main() {
    Fluid fluid(1);  // 1 == Lid-driven cavity
    fluid.setup();
    fluid.run();
    fluid.write_results();

    return 0;
}