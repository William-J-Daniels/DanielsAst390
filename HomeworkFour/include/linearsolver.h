#ifndef LA_LINEARSOLVER_H
#define LA_LINEARSOLVER_H

#include <matrix.h>
#include <vector>
#include <cassert>
#include <thread>
#include <algorithm>
#include <cmath>
#include <limits>

#include <iostream>

namespace la {

template <class T>
class LinearSolver
{
    // solves linear systems of the form  Ax=b for x
    // takes a reference to the data, if you need the original make a copy
public:
    LinearSolver() = default;
    LinearSolver<T>(Matrix<T> init_A, std::vector<T> init_b);

    std::vector<T> gaus_elim();
    std::vector<T> jacobi_iter(std::vector<T> InitGuess);

private:
    Matrix<T> A;
    std::vector<T> b;
    std::vector<T> Solutions;

    unsigned num_pivots = 0;

    std::vector<T> find_scales();
    void find_scales_helper(unsigned start, unsigned end,
                            std::vector<T>& Scales);
    void swap_rows(std::size_t north, std::size_t south);
    void swap_rows(std::size_t north, std::size_t south,
                   std::vector<T> extra_vec);
    void forward_eliminate(std::size_t startr, std::size_t endr,
                           std::size_t current_row);
    void back_substitution();
};

template <class T>
LinearSolver<T>::LinearSolver(Matrix<T> init_A, std::vector<T> init_b)
{
    assert(init_A.numrows() == init_b.size());
    assert(init_A.numrows() == init_A.numcolumns()); // only square for now

    A = init_A;
    b = init_b;
}

template <class T>
std::vector<T> LinearSolver<T>::gaus_elim()
{
    Solutions = std::vector<T>(A.numcolumns());
    auto Scales = find_scales();
    double scaled_col, col_max, row_max_idx;
    std::size_t mod, start, end; // for distributing work to threads

    for (std::size_t current_row = 0; current_row < A.numrows(); current_row++)
    { // loop over all rows
        std::vector<std::thread> Threads;
        // first we loop for pivoting opperturnities
        col_max = std::numeric_limits<T>::lowest();
        row_max_idx = current_row;
        for (std::size_t prow = current_row; prow < A.numrows(); prow++)
        {
            scaled_col = std::abs(A(prow, current_row)) / Scales[prow];
            if (scaled_col > col_max)
            {
                col_max = scaled_col;
                row_max_idx = prow;
            }
        }
        if (row_max_idx-current_row)
            swap_rows(row_max_idx, current_row, Scales);

        // then we forward eliminate
        // forward_eliminate(0, A.numcolumns(), current_row);

        /* parallelize later */
        for (std::size_t erow = current_row + 1; erow < A.numrows(); ++erow)
        {
            auto factor = A(erow, current_row) / A(current_row, current_row);
            for (std::size_t ecol = current_row + 1; ecol < A.numcolumns(); ++ecol)
            {
                A.set(erow, ecol, A(erow,ecol) - A(current_row, ecol) * factor);
            }
            A.set(erow, current_row, 0.0);
            b[erow] = b[erow] - b[current_row]*factor;
        }

        // this bit of code apears multiple times across this class and matrix,
        // it should be moved to a function
        // for (std::size_t t = 0; t < std::thread::hardware_concurrency(); t++)
        // {
        //     // load balance is less than ideal, but better approaches are too
        //     // sophisticated for me
        //     // https://cs.wmich.edu/elise/courses/cs626/s12/cs6260_presentation_1.pdf
        //     mod = (A.numrows()-current_row)%std::thread::hardware_concurrency();
        //     start = (A.numrows()-current_row)/std::thread::hardware_concurrency()*t;
        //     end = (A.numrows()-current_row)/std::thread::hardware_concurrency()*(t+1);
        //     if(mod)
        //     {
        //         if (t < mod)
        //         {
        //             start += t;
        //             end   += t+1;
        //         } else {
        //             start += mod;
        //             end   += mod;
        //         }
        //     }
        //     if (start != end)
        //     {
        //         std::cout << "make thread on row" << current_row << std::endl;
        //         Threads.push_back(std::thread(
        //             &LinearSolver::forward_eliminate, this, start, end,
        //             current_row
        //         ));
        //     }
        // }
        // for (auto & tr : Threads)
        //     tr.join();

    }
    back_substitution();

    return Solutions;
}

template <class T>
std::vector<T> LinearSolver<T>::find_scales()
{
    // returns a vector containing the largest magnitude of each row
    std::vector<T> Scales(A.numrows());
    std::vector<std::thread> Threads;

    auto mod = A.numrows()%std::thread::hardware_concurrency();
    std::size_t start, end;
    for (std::size_t i = 0; i < std::thread::hardware_concurrency(); i++)
    {
        start = (A.numrows()/std::thread::hardware_concurrency())*i;
        end   = (A.numrows()/std::thread::hardware_concurrency())*(i+1);
        if(mod)
        {
            if (i < mod)
            {
                start += i;
                end   += i+1;
            } else {
                start += mod;
                end   += mod;
            }
        } else {
            start += mod*(i+1);
            end   += mod*(i+1);
        }
        if (start != end)
        {
            Threads.push_back(std::thread(
                &LinearSolver::find_scales_helper, this, start, end,
                std::ref(Scales)
            ));
        }
    }
    for (auto& t : Threads)
        t.join();

    return Scales;
}

template <class T>
void LinearSolver<T>::find_scales_helper(unsigned start, unsigned end,
                                         std::vector<T>& Scales)
{
    for (unsigned i = start; i < end; i++)
    {
        Scales[i] = *std::max_element(
            A.slice_begin(i,0,i,A.numcolumns()-1),
            A.slice_end(i,0,i,A.numcolumns()-1),
            [] (T a, T b) { return std::abs(a) < std::abs(b); }
        );
    }
}

template <class T>
void LinearSolver<T>::swap_rows(std::size_t north, std::size_t south)
{
    // can be parallelized later
    std::vector<T> RowTemp(A.numcolumns());
    T VecTemp;

    std::copy(
        A.slice_begin(north,0,north,A.numcolumns()),
        A.slice_end(north,0,north,A.numcolumns()),
        RowTemp.begin()
    );

    std::copy(
        A.slice_begin(south,0,south,A.numcolumns()),
        A.slice_end(south,0,south,A.numcolumns()),
        A.slice_begin(north,0,north,A.numcolumns())
    );

    std::copy(
        RowTemp.begin(),
        RowTemp.end(),
        A.slice_begin(south,0,south,A.numcolumns())
    );

    VecTemp = b[north];
    b[north] = b[south];
    b[south] = VecTemp;

    num_pivots++;
}

template <class T>
void LinearSolver<T>::swap_rows(std::size_t north, std::size_t south,
                   std::vector<T> extra_vec)
{
    T VecTemp;

    swap_rows(north, south);

    VecTemp = extra_vec[north];
    extra_vec[north] = extra_vec[south];
    extra_vec[south] = VecTemp;
}

template <class T>
void LinearSolver<T>::forward_eliminate(std::size_t startr, std::size_t endr,
                                        std::size_t current_row)
{
    T factor;
    for (std::size_t i = startr; i < endr; i++)
    {
        factor = A(i, current_row) / A(current_row, current_row);
        //std::cout << " " << startr << " " << factor << std::endl;
        // std::cout << factor << std::endl;
        std::transform(
            A.slice_begin(i, i, i, A.numcolumns()-1),
            A.slice_end(i, i, i, A.numcolumns()-1),
            A.slice_begin(i, i, i, A.numcolumns()-1),
            [&] (T x) { return x - factor; }
        );
        /* almost there, just need to get the right factor and switch the
         * transform to subtract. need to use two slices on transform
         */
    }
}

template <class T>
void LinearSolver<T>::back_substitution()
{
    Solutions[A.numcolumns()-1] = b[A.numcolumns()-1]/A(A.numcolumns()-1, A.numcolumns()-1);

    for (int brow = (A.numcolumns()-2); brow >= 0; --brow)
    {
        double sum = b[brow];
        for (std::size_t bcol = brow+1; bcol < A.numcolumns(); ++bcol)
        {
            sum = sum - A(brow, bcol) * Solutions[bcol];
        }
        Solutions[brow] = sum / A(brow, brow);
    }
}

template <class T>
std::vector<T> LinearSolver<T>::jacobi_iter(std::vector<T> InitGuess)
{
    assert(A.numrows() == A.numcolumns());
    assert(InitGuess.size() == A.numrows());

    // would be nice to know how to treat sparse matrices here
    Matrix<T> D(A.numrows(), A.numcolumns());
    Matrix<T> R(A.numrows(), A.numcolumns());
    Matrix<T> B(A.numrows(), 1);
    Matrix<T> X(A.numrows(), 1);
    Matrix<T> lastX(A.numrows(), 1);
    for (std::size_t row = 0; row < A.numrows(); row++)
    {
        for (std::size_t col = 0; col < A.numcolumns(); col++)
        {
            if (row == col)
                D.set(row, col, 1.0/A(row, col));
            else
                R.set(row, col, A(row, col));
        }
        B.set(row, 0, b[row]);
        X.set(row, 0, InitGuess[row]);
    }

    do
    {
        lastX = X;
        X = R * X;
        X = B - X;
        X = D * X;
        // std::cout << X << std::endl;
    } while (true);

    return Solutions;
}

} // namespace la

#endif // LA_LINEARSOLVER_H
