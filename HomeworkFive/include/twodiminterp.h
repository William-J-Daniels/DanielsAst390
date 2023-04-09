#ifndef WILL_TWODIMINTERP_H
#define WILL_TWODIMINTERP_H

#include <vector>
#include <functional>
#include <cassert>
#include <cmath>

namespace will {

class TwoDimInterp
{
public:
    TwoDimInterp(std::function<double(double, double)> init_Func,
                 unsigned rows, double start_row, double end_row,
                 unsigned cols, double start_col, double end_col,
                 double base);

    double value_at(double rowval, double colval);

private:
    std::vector<double> Data, row_inputs, col_inputs;
    std::function<double(double, double)> Func;
};

} // namespace will

#endif // WILL_TWODIMINTERP_H
