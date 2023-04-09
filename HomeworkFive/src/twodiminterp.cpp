#include <twodiminterp.h>
#include <logspace.h>

using namespace will;

TwoDimInterp::TwoDimInterp(std::function<double(double, double)> init_Func,
                           unsigned rows, double start_row, double end_row,
                           unsigned cols, double start_col, double end_col,
                           double base)
: Func(init_Func)
{
    assert(rows > 0 && cols > 0);
    assert(start_row < end_row && start_col < end_col);
    assert(base > 1.0);

    // initialize storage
    Data = std::vector<double> (rows*cols);
    row_inputs = std::vector<double> (rows);
    col_inputs = std::vector<double> (cols);

    // fill vectors with input values-- not easily parallelizable, resembles a
    // markov chain
    std::generate_n(
        row_inputs.begin(),
        rows,
        will::Logspace<double> (
            std::pow(base, start_row),
            std::pow(base, (end_row - start_row) / rows)
        )
    );
    std::generate_n(
        col_inputs.begin(),
        cols,
        will::Logspace<double> (
            std::pow(base, start_col),
            std::pow(base, (end_col - start_col) / rows)
        )
    );

    // fill the table-- embarrassingly parallelizable, use omp
    #pragma omp parallel for
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            Data[i*cols + j] = Func(row_inputs[i], col_inputs[j]);
        }
    }
}

double TwoDimInterp::value_at(double rowval, double colval)
{
    assert(rowval > row_inputs.front() && rowval < row_inputs.back());
    assert(colval > col_inputs.front() && colval < col_inputs.back());

    // find grid points surrounding point of interest
    auto lower_row = std::lower_bound(row_inputs.begin(), row_inputs.end(),
                                      rowval);
    assert(std::distance(lower_row, row_inputs.end()) > 1);
    double lrow = *lower_row, hrow = *(++lower_row);

    auto lower_col = std::lower_bound(col_inputs.begin(), col_inputs.end(),
                                      colval);
    assert(std::distance(lower_row, row_inputs.end()) > 1);
    double lcol = *lower_col, hcol = *(++lower_col);

    // compute coefficients
    double a = (
        Data[((lower_row+1) - row_inputs.begin()) * col_inputs.size() +
             ((lower_col+1) - col_inputs.begin())] -
        Data[(lower_row - row_inputs.begin()) * col_inputs.size() +
             (lower_col - col_inputs.begin())]
    ) / (
        std::sqrt(std::pow(hrow-lrow, 2) + std::pow(hcol-lcol, 2))
    );

    double b = (
        Data[((lower_row+1) - row_inputs.begin()) * col_inputs.size() +
             (lower_col - col_inputs.begin())] -
        Data[((lower_row) - row_inputs.begin())   * col_inputs.size() +
             (lower_col - col_inputs.begin())]
    ) / (
        hrow - lrow
    );

    double c = (
        Data[(lower_row - row_inputs.begin()) * col_inputs.size() +
             ((lower_col+1) - col_inputs.begin())] -
        Data[(lower_row - row_inputs.begin()) * col_inputs.size() +
             (lower_col - col_inputs.begin())]
    ) / (
        hcol - lcol
    );

    double d = Data[(lower_row - row_inputs.begin()) * col_inputs.size() +
                    (lower_col - col_inputs.begin())];

    // return value at point of interest
    return(
        a * (
            rowval - lrow
        ) * (
            colval - lcol
        ) + b * (
            rowval - lrow
        ) + c * (
            colval - lcol
        ) + d
    );
}
