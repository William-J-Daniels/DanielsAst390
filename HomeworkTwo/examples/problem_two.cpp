#include <simpsons.h>
#include <indefinite_integrals.h>
#include <cmath>
#include <iostream>
#include <limits>

// constants
const double kB            = 1.380649e-23;
const double T             = 1.5e7;
const double mp            = 1.67262492369e-27;
const double analtic       = 2.0 * std::sqrt((2.0*kB*T)/(M_PI*mp));
const double factor        = 4.0 * std::sqrt((2.0 * kB * T) / (mp * M_PI));
const double integrand_max = 1.225;

// function prototypes
double integrand(double input);


int main()
{
    auto Integral = will::Simpsons(integrand);

    double result = Integral.converge(0.0, 1.0, 2,
                                      std::numeric_limits<unsigned>::max());
    result *= factor * integrand_max;

    std::cout << "Analytic: " << analtic << std::endl
              << "Numerical: " << result << std::endl
              << "Iterations: " << Integral.get_num_iter() << std::endl
              << "Difference : " << std::abs(result - analtic) << std::endl;


    return 0;
}


// function definitions
double integrand(double input)
{
    double transformed_input = will::infinite_transform(input, integrand_max);
    return(
        std::exp(-std::pow(transformed_input, 2)) *
        std::pow(transformed_input, 3) *
        std::pow(1.0 - input + will::SMALL, -2)
    );
}
