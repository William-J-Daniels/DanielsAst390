#include <twodiminterp.h>
#include <logspace.h>
#include <iostream>
#include <cmath>
#include <iostream>
#include <vector>

double pressure(double density, double temperature);

int main()
{
    for (unsigned i = 6; i < 15; i++) // not enough memory to go bigger
    {
        unsigned size = std::pow(2, i);
        std::cout << size << " by " << size << " grid" << std::endl;
        auto PTab = will::TwoDimInterp(pressure,
                                        size, 0.0, 4.0,
                                        size, 6.0, 8.0,
                                        10.0);
        auto interp = PTab.value_at(500.0, 6.0e7);
        auto anal   = pressure(500.0, 6.0e7);
        std::cout << "Interpolated: " << interp << std::endl;
        std::cout << "Analytic: " << anal << std::endl;
        std::cout << "Difference: " << std::abs(interp - anal) << std::endl
                  << std::endl;
    }

    return 0;
}

double pressure(double density, double temperature)
{
    assert(density > 0 && temperature > 0);
    double kb = 1.38e-16, // K
           a  = 7.56e-15, // erg cm^-3
           ui = 1.26,     // dimless
           ue = 1.15,     // dimless
           k  = 10.0e13;  // erg cm^-3

    return(
        (
            kb * density * temperature
        ) / (
            ui * ue
        ) * (
            a * std::pow(temperature, 4)
        ) / (
            3.0
        ) + (
            k *
            std::pow(density, 5.0/3.0) *
            std::pow(ue, -5.0/3.0)
        )
    );
}
