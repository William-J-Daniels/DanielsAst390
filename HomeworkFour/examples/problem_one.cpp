#include <matrix.h>
#include <linearsolver.h>
#include <iostream>
#include <vector>

int main()
{
    for (unsigned i = 2; i <= 15; i++)
    {
        auto H = la::Matrix<double>(i, i);
        auto x = std::vector<double>(i);

        for (unsigned j = 0; j < i; j++)
        {
            for (unsigned k = 0; k < i; k++)
            {
                H.set(j,k, 1.0/(j+k+1.0));
            }
            x[j] = static_cast<double>(j);
        }

        auto b = H * x;
        auto Solver = la::LinearSolver(H, b);
        auto Solution = Solver.gaus_elim();
        auto Errors = std::vector<double>(i);
        std::transform(
            Solution.begin(), Solution.end(),
            x.begin(), Errors.begin(),
            [&](double a, double b){return std::abs(b - a);}
        );
        std::cout << "Size of H: " << i << "; Maximum error: " <<
                     *std::max_element(Errors.begin(), Errors.end()) <<
                     std::endl;
    }

    return 0;
}
