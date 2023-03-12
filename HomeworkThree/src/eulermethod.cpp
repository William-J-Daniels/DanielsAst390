#include <eulermethod.h>


#include <iostream>

using namespace will;

std::pair<double, double> EulerMethod::advance(double step)
{
    std::get<1>(state) += step * DiffEq(std::get<1>(state));
    std::get<0>(state) += step;

    return state;
}
