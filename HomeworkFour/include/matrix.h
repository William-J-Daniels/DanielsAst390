#ifndef LA_MATRIX
#define LA_MATRIX

#include <vector>
#include <cassert>
#include <iterator>
#include <cstddef>
#include <algorithm>
#include <iostream>
#include <thread>
#include <cstdlib>
#include <execution>

namespace la {

template <class T>
class Matrix
{
public:
    // constructors
    Matrix() = default;
    Matrix(std::size_t init_rows, std::size_t init_columns);
    Matrix(std::size_t init_rows, std::size_t init_columns, T init_val);
    Matrix(std::vector<std::vector<T>> v);

    // accessors and mutators
    std::size_t numrows() { return rows; }
    std::size_t numcolumns() { return columns; }
    T operator() (std::size_t row, std::size_t column);
    friend std::ostream& operator<< (std::ostream& out, const Matrix<T>& M)
    {
        for (int row = 0; row < M.rows; row++)
        {
            for (int column = 0; column < M.columns; column++)
            {
                out << M.data[row*M.columns + column] << " ";
            }
            out << std::endl;
        }
        return out;
    }

    void set(std::size_t row, std::size_t column, T val);

    // iterators
    // row iterator is just the iterator of the internal vector
    auto row_begin() {return data.begin();}
    auto row_end()   {return data.end();}
    class columnIterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = value_type*;
        using reference         = value_type&;

        columnIterator(pointer ptr, Matrix &Owner) : m_ptr(ptr), M(Owner) {}
        reference operator*() const {return *m_ptr;}
        pointer operator->() {return m_ptr;}
        columnIterator& operator++()
        {
            if (counter%M.columns != 0)
            { // when not at the end
                // stride by one row
                m_ptr = m_ptr + M.rows;
                counter = counter + 1;
                return *this;
            }
            if (counter == M.rows*M.columns)
            { // if we've visited all the addresses
                // stride by one address
                m_ptr = m_ptr + 1;
                return *this;
            }
            // when at the end of a column
            // stride to the next column
            m_ptr = m_ptr - M.rows * (M.columns - 1) + 1;
            counter = counter + 1;
            return *this;
        }
        columnIterator operator++(int)
            {columnIterator tmp = *this; ++(*this); return tmp;}
        friend bool operator== (const columnIterator& a,const columnIterator& b)
            {return a.m_ptr == b.m_ptr;}
        friend bool operator!= (const columnIterator& a,const columnIterator& b)
            {return a.m_ptr != b.m_ptr;}
    private:
        pointer m_ptr;
        Matrix<T> M;
        unsigned counter = 1;
    };
    auto column_begin() {return columnIterator(&data[0],            *this);}
    auto column_end()   {return columnIterator(&data[rows*columns], *this);}
    class Slice
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = value_type*;
        using reference         = value_type&;

        Slice(pointer ptr, Matrix &Owner,
              std::size_t start_row, std::size_t start_col,
              std::size_t end_row, std::size_t end_col) :
              m_ptr(ptr), M(Owner),
              sr(start_row), sc(start_col),
              er(end_row), ec(end_col) {}
        reference operator*() const {return *m_ptr;}
        pointer operator->() {return m_ptr;}
        Slice& operator++()
        {
            if (counter%((ec-sc)+1) == 0)
            { // when we're not the end of the slice's row
                // get to the next row of the slice
                m_ptr += M.columns-(ec-sc);
                counter += 1;
                return *this;
            }
            // default behavior is to advance by one
            m_ptr += 1;
            counter += 1;
            return *this;
        }
        Slice operator++(int)
            {Slice tmp = *this; ++(*this); return tmp;}
        friend bool operator== (const Slice& a,const Slice& b)
            {return a.m_ptr == b.m_ptr;}
        friend bool operator!= (const Slice& a,const Slice& b)
            {return a.m_ptr != b.m_ptr;}
    private:
        pointer m_ptr;
        Matrix<T> M;
        size_t sr, sc, er, ec;
        unsigned counter = 1;
    };
    auto slice_begin(std::size_t start_row, std::size_t start_col,
                     std::size_t end_row,   std::size_t end_col)
    {
        return Slice(&data[start_row*rows+start_col], *this,
                     start_row, start_col, end_row, end_col);
    }
    auto slice_end(std::size_t start_row, std::size_t start_col,
                   std::size_t end_row,   std::size_t end_col)
    {
        return Slice(
            &data[(end_row+1)*rows+(end_col-(end_col-start_col))], *this,
            start_row, start_col, end_row, end_col);
    }

    // matrix operations (will grow over time)
    Matrix<T> operator- ();
    Matrix<T> operator+ (Matrix<T>& M2);
    Matrix<T> operator- (Matrix<T>& M2);

    friend Matrix<T> operator* (T scalar, Matrix<T>& M)
    {
        auto newMatrix = Matrix(M.rows, M.columns);

        std::transform(
            std::execution::par,
            M.row_begin(), M.row_end(),
            newMatrix.row_begin(),
            [&scalar] (T a) { return a * scalar; }
        );

        return newMatrix;
    }
    friend Matrix<T> operator* (Matrix<T>& M, T scalar)
    {
        auto newMatrix = Matrix(M.rows, M.columns);

        std::transform(
            std::execution::par,
            M.row_begin(), M.row_end(),
            newMatrix.row_begin(),
            [&scalar] (T a) { return a * scalar; }
        );

        return newMatrix;
    }
    void transpose();
    Matrix<T> operator* (Matrix<T>& M2);
    Matrix<T> naive_mult(Matrix<T>& M2);
    Matrix<T> strassen_mult(Matrix<T>& M2);

private:
    std::size_t rows, columns;
    std::vector<T> data;
    bool row_major = true; // always stored row major, really refers to whether
                           // the cache is setup for multiplication

    void mult_helper(Matrix<T>& M2, Matrix<T>& Mout,
                     size_t thread_idx);
};

template<class T>
Matrix<T>::Matrix(std::size_t init_rows, std::size_t init_columns)
{
    assert(init_rows > 0 && init_columns > 0);

    rows = init_rows;
    columns = init_columns;

    data = std::vector<T>(rows*columns);
}

template <class T>
Matrix<T>::Matrix(std::size_t init_rows, std::size_t init_columns, T init_val)
{
    assert(init_rows > 0 && init_columns > 0);

    rows = init_rows;
    columns = init_columns;

    data = std::vector<T>(rows*columns, init_val);
}

template <class T>
Matrix<T>::Matrix(std::vector<std::vector<T>> v)
{ // https://zingale.github.io/computational_astrophysics/basics/linear-algebra/la-cxx.html
    rows    = v.size();
    columns = v[0].size();
    data.resize(rows*columns);

    std::size_t idx = 0;

    for (std::size_t i = 0; i < rows; ++i)
    {
        assert(v[i].size() == columns);
        for (std::size_t j = 0; j < columns; ++j)
        {
            data[idx] = v[i][j];
            ++idx;
        }
    }
}

template <class T>
T Matrix<T>::operator() (std::size_t row, std::size_t column)
{
    assert(row >= 0 && row < rows);
    assert(column >= 0 && column < columns);

    return data[row*columns + column];
}

template <class T>
void Matrix<T>::set(std::size_t row, std::size_t column, T val)
{
    data[row*columns + column] = val;
}

// matrix operations

template <class T>
Matrix<T> Matrix<T>::operator-()
{
    auto newMatrix = Matrix<T>(rows, columns);
    std::transform(
        std::execution::par,
        row_begin(), row_end(),
        newMatrix.row_begin(),
        std::negate<>{}
    );

    return newMatrix;
}

template <class T>
Matrix<T> Matrix<T>::operator+ (Matrix<T>& M2)
{
    // only defined for matrices of the same datatype
    assert(rows == M2.rows && columns == M2.columns);

    auto newMatrix = Matrix<T>(rows, columns);
    std::transform(
        std::execution::par,
        row_begin(), row_end(),
        M2.row_begin(),
        newMatrix.row_begin(),
        std::plus<>{}
    );

    return newMatrix;
}

template <class T>
Matrix<T> Matrix<T>::operator- (Matrix<T>& M2)
{
    assert(rows == M2.rows && columns == M2.columns);

    auto newMatrix = Matrix<T>(rows, columns);
    std::transform(
        std::execution::par,
        row_begin(), row_end(),
        M2.row_begin(),
        newMatrix.row_begin(),
        std::minus<>{});

    return newMatrix;
}

template <class T>
void Matrix<T>::transpose()
{ // https://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c
    // this method changes the data- change is noted with row_major variable
    // can imporve parallelization later

    auto TempData = std::vector<T>(data.size());

    std::copy(
        std::execution::par,
        data.begin(), data.end(),
        TempData.begin()
    );

    for (int n = 0; n < rows*columns; n++)
    {
        auto dv = std::div(n, rows);

        data[n] = TempData[columns*dv.rem + dv.quot];
    }

    std::swap(rows, columns);
    row_major = false;
}

template <class T>
Matrix<T> Matrix<T>::operator*(Matrix<T>& M2)
{
    assert(columns = M2.rows);

    auto newMatrix = Matrix<T>(rows, M2.columns, 0.0);
    /*
     * Huang, Jianyu; Smith, Tyler M.; Henry, Greg M.; van de Geijn, Robert A.
     * (13 Nov 2016) found that Strassan's algorithm becomes useful at about
     * 512x512.
     */

    if (rows < 511)
    {
        newMatrix = naive_mult(M2);
    } else {
        newMatrix = strassen_mult(M2);
    }
return newMatrix;
}

template <class T>
Matrix<T> Matrix<T>::naive_mult(Matrix<T>& M2)
{
    std::vector<std::thread> Threads;
    auto newMatrix = Matrix<T>(rows, M2.columns);

    // if (M2.rows == M2.columns && M2.row_major)
    //     M2.transpose(); // figure out later have HW to do

    for (unsigned i = 0; i < std::thread::hardware_concurrency(); i++)
    {
        Threads.push_back(
            std::thread(&Matrix<T>::mult_helper, this,
                        std::ref(M2), std::ref(newMatrix), i)
        );
    }

    for (auto& t : Threads)
        t.join();

    return newMatrix;
}

template <class T>
void Matrix<T>::mult_helper(Matrix<T>& M2, Matrix<T>& Mout,
                            size_t thread_idx)
{
    double sum;
    auto remainder = rows%std::thread::hardware_concurrency();

    std::size_t start = (rows/std::thread::hardware_concurrency())*thread_idx;
    std::size_t end = (rows/std::thread::hardware_concurrency())*(thread_idx+1);
    if(remainder)
    {
        // when the number of rows is not divisible by the number of threads...
        if (thread_idx < remainder)
        {
            // and the thread needs to do an extra loop...
            // add the index to the start, because this is equal to the number
            // of threads before it that had to do an extra loop,
            start += thread_idx;
            // and increase the end by one plus the index to give it more work.
            end   += thread_idx+1;
        } else {
            // and the thread does not need to do an extra loop...
            // add the index to the start, because this is equal to the number
            // of threads before it that had to do an extra loop,
            start += thread_idx;
            // and add the index to the end, becuase this thread doesn't need to
            // do extra work.
            end   += thread_idx;
        }
    } else {
        // when the number of rows is divisible by the number of threads, divide
        // the work trivially
        start += remainder*(thread_idx+1);
        end   += remainder*(thread_idx+1);
    }


    for (std::size_t i = start; i < end; i++)
    {
        for (std::size_t j = 0; j < M2.columns; j++)
        {
            sum = 0.0;
            for (std::size_t k = 0; k < columns; k++)
            {
                sum += data[i*columns+k] * M2(k,j);
            }
            Mout.set(i,j, sum);
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::strassen_mult(Matrix<T>& M2)
{
    /*
     * An obvious way to thread Strassen's method is to create a thread for each
     * of the 7 child multiplications. On my system, this is unideal becuase I
     * only have two physical cores and four logcal threads. The result of this
     * would be multiple naive multiplications running in parallel, but this is
     * no more efficient than performing each multiplication with a threaded
     * implimentation serially, as long as I take advantage of std::execution to
     * make the copying and addition parallel. Might incur additional thread
     * overhead, but for large matrices this won't matter and this approach is
     * best for portability amoung systems with different numbers of CPUs
     */
    assert(rows == columns && rows%2 == 0 && rows == columns);
    // will impliment padding later
    auto newMatrix = Matrix<T>(rows, M2.columns);

    std::size_t size = rows/2;
    std::vector<Matrix<T>> Intermediates(14, Matrix(size, size));
    std::vector<Matrix<T>> Products(7, Matrix(size, size));

    // make all the intermediate sums
    // can be parallelized-- maybe later, if we make a thread pool for class
    std::transform( // A11 + A22
        std::execution::par,
        slice_begin(0,0,size-1,size-1), slice_end(0,0,size-1,size-1),
        slice_begin(size, size, rows-1, columns-1), Intermediates[0].row_begin(),
        std::plus<>{}
    );
    std::transform( // B11 + B22
        std::execution::par,
        M2.slice_begin(0,0,size-1,size-1), M2.slice_end(0,0,size-1,size-1),
        M2.slice_begin(size, size, rows-1, columns-1), Intermediates[1].row_begin(),
        std::plus<>{}
    );
    std::transform( // A21 + A22
        std::execution::par,
        slice_begin(size,0,rows-1,size-1), slice_end(size,0,rows-1,size-1),
        slice_begin(size, size, rows-1, columns-1), Intermediates[2].row_begin(),
        std::plus<>{}
    );
    std::transform( // B12 - B22
        std::execution::par,
        M2.slice_begin(0,size,size-1,columns-1), M2.slice_end(0,size,size-1,columns-1),
        M2.slice_begin(size, size, rows-1, columns-1), Intermediates[3].row_begin(),
        std::minus<>{}
    );
    std::transform( // B21 - B11
        std::execution::par,
        M2.slice_begin(size,0,columns-1,size-1), M2.slice_end(size,0,columns-1,size-1),
        M2.slice_begin(0,0,size-1,size-1), Intermediates[4].row_begin(),
        std::minus<>{}
    );
    std::transform( // A11 + A12
        std::execution::par,
        slice_begin(0,0,size-1,size-1), slice_end(0,0,size-1,size-1),
        slice_begin(0,size,size-1,rows-1), Intermediates[5].row_begin(),
        std::plus<>{}
    );
    std::transform( // A21 - A11
        std::execution::par,
        slice_begin(size,0,columns-1,size-1), slice_end(size,0,columns-1,size-1),
        slice_begin(0,0,size-1,size-1), Intermediates[6].row_begin(),
        std::minus<>{}
    );
    std::transform( // B11 + B12
        std::execution::par,
        M2.slice_begin(0,0,size-1,size-1), M2.slice_end(0,0,size-1,size-1),
        M2.slice_begin(0,size,size-1,rows-1), Intermediates[7].row_begin(),
        std::plus<>{}
    );
    std::transform( // A12 - A22
        std::execution::par,
        slice_begin(0,size,size-1,rows-1), slice_end(0,size,size-1,rows-1),
        slice_begin(size,size,rows-1,columns-1), Intermediates[8].row_begin(),
        std::minus<>{}
    );
    std::transform( // B21 + B22
        std::execution::par,
        M2.slice_begin(size,0,columns-1,size-1), M2.slice_end(size,0,columns-1,size-1),
        M2.slice_begin(size,size,rows-1,columns-1), Intermediates[9].row_begin(),
        std::plus<>{}
    );
    std::copy( // B11
        std::execution::par,
        M2.slice_begin(0,0,size-1,size-1), M2.slice_end(0,0,size-1,size-1),
        Intermediates[10].row_begin()
    );
    std::copy( // A11
        std::execution::par,
        slice_begin(0,0,size-1,size-1), slice_end(0,0,size-1,size-1),
        Intermediates[11].row_begin()
    );
    std::copy( // A22
        std::execution::par,
        slice_begin(size,size,rows-1,columns-1), slice_end(size,size,rows-1,columns-1),
        Intermediates[12].row_begin()
    );
    std::copy( // B22
        std::execution::par,
        M2.slice_begin(size,size,rows-1,columns-1), M2.slice_end(size,size,rows-1,columns-1),
        Intermediates[13].row_begin()
    );

    // perform multiplications don't parallelize-- naive_mult already is
    Products[0] = Intermediates[0] * Intermediates[1];
    Products[1] = Intermediates[2] * Intermediates[10];
    Products[2] = Intermediates[11] * Intermediates[3];
    Products[3] = Intermediates[12] * Intermediates[4];
    Products[4] = Intermediates[5] * Intermediates[13];
    Products[5] = Intermediates[6] * Intermediates[7];
    Products[6] = Intermediates[8] * Intermediates[9];

    // transform products into the new matrix
    // can also be put into a thread pool
    std::transform( // C11 = M1 + M4
        std::execution::par,
        Products[0].row_begin(), Products[0].row_end(),
        Products[3].row_begin(), newMatrix.slice_begin(0,0,size-1,size-1),
        std::plus<>{}
    );
    std::transform( // C11 = C11 - M5 (C11 = M1 + M4 - M5)
        std::execution::par,
        newMatrix.slice_begin(0,0,size-1,size-1), newMatrix.slice_end(0,0,size-1,size-1),
        Products[4].row_begin(), newMatrix.slice_begin(0,0,size-1,size-1),
        std::minus<>{}
    );
    std::transform( // C11 = C11 + M7 (C11 = M1 + M4 - M5 + M7)
        std::execution::par,
        newMatrix.slice_begin(0,0,size-1,size-1), newMatrix.slice_end(0,0,size-1,size-1),
        Products[6].row_begin(), newMatrix.slice_begin(0,0,size-1,size-1),
        std::plus<>{}
    );
    std::transform( // C12 = M3 + M5
        std::execution::par,
        Products[2].row_begin(), Products[2].row_end(),
        Products[4].row_begin(), newMatrix.slice_begin(0,size,size-1,columns-1),
        std::plus<>{}
    );
    std::transform( // C21 = M2 + M4
        std::execution::par,
        Products[1].row_begin(), Products[1].row_end(),
        Products[3].row_begin(), newMatrix.slice_begin(size,0,rows-1,size-1),
        std::plus<>{}
    );
    std::transform( // C22 = M1 - M2
        std::execution::par,
        Products[0].row_begin(), Products[0].row_end(),
        Products[1].row_begin(), newMatrix.slice_begin(size,size,rows-1,columns-1),
        std::minus<>{}
    );
    std::transform( // C22 = C22 + M3 (C22 = M1 - M2 + M3)
        std::execution::par,
        newMatrix.slice_begin(size,size,rows-1,columns-1), newMatrix.slice_end(size,size,rows-1,columns-1),
        Products[2].row_begin(), newMatrix.slice_begin(size,size,rows-1,columns-1),
        std::plus<>{}
    );
    std::transform( // C22 = C22 + M6 (C22 = M1 - M2 + M3 + M6)
        std::execution::par,
        newMatrix.slice_begin(size,size,rows-1,columns-1), newMatrix.slice_end(size,size,rows-1,columns-1),
        Products[5].row_begin(), newMatrix.slice_begin(size,size,rows-1,columns-1),
        std::plus<>{}
    );

    return newMatrix;
}

} // namespace la

#endif // LA_MATRIX
