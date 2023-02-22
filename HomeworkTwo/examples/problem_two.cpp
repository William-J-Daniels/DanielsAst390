#include <simpsons.h>
#include <cmath>
#include <iostream>
#include <limits>

const double SMALL   = 1.0e-30;
const double kB      = 8.617333262e-5;
const double T       = 1.5e7;
const double mp      = 938.27208816e6;
const double analtic = 2.0 * std::sqrt((2.0*kB*T)/(M_PI*mp));

double zv(double x, double c)
{
    return x / (c + x);
}
double xv(double z, double c)
{
    return c*z / (1.0 - z + SMALL);
}
double integrand(double p)
{
    p = zv(p, 3.0);
    return(
        std::pow(2.0 * M_PI * mp * kB * T, -3.0/2.0) *
        std::exp(-std::pow(p, 2) / (2 * mp * kB * T)) *
        4.0 * M_PI * std::pow(p, 2) *
        (p / mp)
    );
}

int main()
{
    auto test = will::Simpsons(integrand);

    std::cout << analtic << std::endl;

    std::cout << test.converge(0.0, 1.0, 1024, std::numeric_limits<unsigned>::max()) << ", " << test.get_num_iter();

    // for (double x = 0.0; x<1.0; x+=0.001)
    //     std::cout<< integrand(x) << std::endl;

    return 0;
}
