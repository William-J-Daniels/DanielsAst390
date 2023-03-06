#include <eulercromer.h>

using namespace will;

std::pair<double, double> EulerCromer::advance(double step)
{
    std::get<0>(state) += step;
    std::get<1>(state) += step * DiffEq(std::get<0>(state));

    return state;
}
