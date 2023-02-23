#ifndef WILL_SIMPSONS_H
#define WILL_SIMPSONS_H

#include <functional>
#include <limits>

namespace will {

class Simpsons
{
public:
    Simpsons(std::function<double(double)> init_Function);

    double evaluate(double start, double end, unsigned intervals);
    double converge(double start, double end, unsigned intervals,
                    unsigned max_iter);

    void reset();
    double get_integral();
    unsigned get_num_iter();

private:
    std::function<double(double)> Function;
    double integral = 0.0;
    unsigned num_iter = 1;
    double precision = 1.0e-8;
};

} // namespace will

#endif // WILL_SIMPSONS_H
