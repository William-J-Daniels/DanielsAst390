#ifndef DIFFERENTIATOR_H
#define DIFFERENTIATOR_H

#include <functional>
#include <iostream>
#include <cassert>
#include <cmath>

namespace will{

class Differentiator
{
    // second order numerical approximation of the derivative
    // provides non iterative and iterative methods
public:
    // can construct defaultly or by passing the function
    Differentiator() = delete;
    Differentiator(std::function<double(double)> init_Function);

    // solvers
    double evaluate(double location);
    double evaluate(double location, double new_step_size);
    double converge(double location);
    double converge(double location, double new_precision);

    // mutators
    void set_step_size(double new_step_size);
    void reset();

private:
    std::function<double(double)> Function; // function to be differentiated
    double step_size = 1.0e-6; // default initail step size
    double precision = 1.0e-12; // precision attempted to converge to
    double result; // the approximated derivative
};

} // namespace will

#endif // DIFFERENTIATOR_H
