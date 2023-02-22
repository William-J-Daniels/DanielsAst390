#include <iostream>
#include <cmath>

// function signatures
double Original(double input);
double Equivalent(double input);

int main()
{
    std::cout << "At 1.0e-6: "
              << Original(1.0e-6) << ", " << Equivalent(1.0e-6)
              << std::endl;
    std::cout << "At 1.0e-7: "
              << Original(1.0e-7) << ", " << Equivalent(1.0e-7)
              << std::endl;
    std::cout << "At 1.0e-8: "
              << Original(1.0e-8) << ", " << Equivalent(1.0e-8)
              << std::endl;

    return 0;
}

// function definitions
double Original(double input)
{
    return (std::sqrt(std::pow(input, 2) + 1.0) + 1.0);
}

double Equivalent(double input)
{
    return (std::pow(input, 2) / std::sqrt(std::pow(input, 2) + 1));
}
