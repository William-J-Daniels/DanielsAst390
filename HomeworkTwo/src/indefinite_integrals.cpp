#include <indefinite_integrals.h>

double will::finite_transform(double input, double max)
{
    return(
        input / (max + input)
    );
}

double will::infinite_transform(double input, double max)
{
    return(
        (max * input) / (1.0 - input + will::SMALL)
    );
}
