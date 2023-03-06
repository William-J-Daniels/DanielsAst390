#include <simpsons.h>
#include <cmath>
#include <iostream>
#include <limits>

double Parabola(double input);

int main()
{
    auto ParabolaIntegral = will::Simpsons(Parabola);

    std::cout << "Parabola evauation" << std::endl;
    for (int i = 1; i < 100; i++)
    {
        std::cout << ParabolaIntegral.evaluate(0.0, 1.0, i) << std::endl;
    }
}

double Parabola(double input)
{
    return std::pow(input, 2);
}
