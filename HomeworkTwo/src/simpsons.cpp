#include <simpsons.h>
#include <iostream>

using namespace will;

Simpsons::Simpsons(std::function<double(double)> init_Function)
{
    Function = init_Function;
}

double Simpsons::evaluate(double start, double end, unsigned intervals)
{
    integral = 0.0;
    double differential = (end - start) / intervals;
    double current;

    // loop as if we always have even intervals
    for (unsigned n = 0; n < intervals; n+=2)
    {
        current = start + n*differential;
        integral += Function(current) +
                    4.0 * Function(current + differential) +
                    Function(current + 2.0*differential);
    }
    integral = integral * differential / 3.0;

    if (intervals%2) // take advantage of int <-> bool, this runs if odd
    { // adjust for overshoot caused by odd intervals
        integral -= differential * (
            5.0 * Function(end) +
            8.0 * Function(end + differential) -
            Function(end + 2.0*differential)
        ) / 12.0;
    }

    return integral;
}

double Simpsons::converge(double start, double end, unsigned int intervals,
                          unsigned max_iter)
{
    double last_integral;
    do
    {
        last_integral = integral;
        evaluate(start, end, intervals*num_iter);
        num_iter++;
    } while (std::abs(last_integral - integral) > precision
             && num_iter < max_iter);

    return integral;
}

void Simpsons::reset()
{
    integral = 0.0;
    num_iter = 1;
}

double Simpsons::get_integral()
{
    return integral;
}

unsigned Simpsons::get_num_iter()
{
    return num_iter;
}
