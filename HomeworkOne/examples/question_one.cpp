#include <differentiator.h>
#include <iostream>
#include <cmath>

// function signatures
double MyFunction(double input);
double AnalyticDerivative(double input);

int main()
{
    will::Differentiator Dsin = will::Differentiator(MyFunction);
    double evaluation_point = 0.1; // radians

    for (int i = 1; i < 11; i++)
    {
        std::cout << Dsin.evaluate(evaluation_point, std::pow(0.5, i)) -
                     AnalyticDerivative(evaluation_point) << std::endl;
    }

    return 0;
}

// function definitions
double MyFunction(double input)
{
    return std::sin(input);
}

double AnalyticDerivative(double input)
{
    return -std::sin(input);
}
