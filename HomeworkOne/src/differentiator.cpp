#include <differentiator.h>

using namespace will;

// constructors

Differentiator::Differentiator(std::function<double(double)> init_Function)
{
    Function = init_Function;
}

// solvers
// using time step

double Differentiator::evaluate(double location)
{
    // returns the derivative approximated non-itteratively using step_size
    return (
        (Function(location + step_size) -
         2.0 * Function(location) +
         Function(location - step_size)
        )
        /
        std::pow(step_size, 2)
    );
}

double Differentiator::evaluate(double location, double new_step_size)
{
    // overload of evaluate that also sets a new step size
    step_size = new_step_size;
    return evaluate(location);
}

// iterative

double will::Differentiator::converge(double location)
{
    // iteratively reduces the step size until precision is achieved
    double last_result;
    result = evaluate(location);
    do
    {
        last_result = result;
        step_size = step_size / 2.0;
        result = evaluate(location);
        std::cout << result << std::endl;
    } while (std::abs(result - last_result) > precision);
    return result;
}

double will::Differentiator::converge(double location, double new_precision)
{
    precision = new_precision;
    return converge(location);
}
// mutators

void Differentiator::set_step_size(double new_step_size)
{
    assert(step_size > 0);
    step_size = new_step_size;
}

void Differentiator::reset()
{
    step_size = 1.0e-6;
    precision = 1.0e-12;
}
