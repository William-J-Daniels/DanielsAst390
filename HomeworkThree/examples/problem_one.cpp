#include <eulerpendulum.h>
#include <eulercromerpendulum.h>
#include <iostream>
#include <list>
#include <array>
#include <fstream>
#include <string>


int main()
{
    const std::pair<double, double> INITIALCONDITIONS = {0.0, (10.0*M_PI/180.0)};
    const double TIME_STEP = 0.5;

    auto EulerPend = will::EulerPendulum(INITIALCONDITIONS,
                                         TIME_STEP);
    auto EulerCromerPend = will::EulerCromerPendulum(INITIALCONDITIONS,
                                                     TIME_STEP);

    std::list<std::array<double, 5>> Results; // linked list to store results
    // time, euler theta, euler nrg, cromer theta, cromer nrg
    const std::string Outfile = "../../HomeworkThree/data/P1.csv";
    // from directory of executable

    for (int i = 0; i < 40; i++)
    {
        auto CurrentE  = EulerPend.advance(TIME_STEP);
        auto currentEC = EulerCromerPend.advance(TIME_STEP);
        // std::cout << EulerPend.get_nrg() << ", " << EulerCromerPend.get_nrg() << std::endl;

        std::cout << std::get<2>(CurrentE) << ", " << std::get<2>(currentEC)
                  << std::endl;

        Results.push_back(std::array<double, 5> {i*TIME_STEP,
                                                 std::get<2>(CurrentE),
                                                 EulerPend.get_nrg(),
                                                 std::get<2>(currentEC),
                                                 EulerCromerPend.get_nrg()});
    }

    return 0;
}
