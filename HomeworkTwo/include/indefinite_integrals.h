#ifndef WILL_INDEFINITEINTEGRALS_H
#define WILL_INDEFINITEINTEGRALS_H

#include <functional>

namespace will{

// constants
constexpr double SMALL = 1.0e-30;

// function prototypes
double finite_transform(double input, double max);
double infinite_transform(double input, double max);

} // namespace will

#endif //WILL_INDEFINITEINTEGRALS_H
