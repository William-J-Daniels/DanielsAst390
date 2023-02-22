#include <simpsons.h>
#include <cmath>
#include <iostream>
#include <limits>

// function prototypes
double Function(double input);

int main()
{
    auto Integrator = will::Simpsons(Function);
    std::vector<unsigned> Intervals = {3, 5, 9, 17, 33};
    const double start = 0.0,
                 end   = 5.0,
                 analytic = -5.0 / (2.0 * M_PI);
    double last_error = std::numeric_limits<double>::infinity(),
           error,
           approximation;

    for (auto i : Intervals)
    {
        approximation = Integrator.evaluate(start, end, i);
        error = std::abs(approximation - analytic);

        std::cout << i << ", "
                  << approximation << ", "
                  << error << ", "
                  << last_error / error << std::endl;

        last_error = error;
    }

    return 0;
}

double Function(double input)
{
    return(
        input * std::sin(2.0 * M_PI * input)
    );
}
