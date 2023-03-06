#include <eulermethod.h>
#include <iostream>

double TestFunc(double x)
{
    return x;
}

int main()
{
    std::pair<double, double> init = {0.0, 0.0};
    auto Test = will::EulerMethod(TestFunc, init, 0.1);

    for (int i = 0; i < 10; i++)
    {
        std::cout << std::get<1>(Test.advance(0.1)) << std::endl;
    }

    return 0;
}
