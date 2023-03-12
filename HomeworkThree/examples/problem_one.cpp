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
    // limited by IO, so prefer constant pesh_back over constant random access
    const std::string FileName = "../../HomeworkThree/data/P1.csv";
    // from directory of executable

    for (int i = 0; i < 40; i++)
    {
        auto CurrentE  = EulerPend.advance(TIME_STEP);
        auto currentEC = EulerCromerPend.advance(TIME_STEP);

        Results.push_back(std::array<double, 5> {i*TIME_STEP,
                                                 std::get<2>(CurrentE),
                                                 EulerPend.get_nrg(),
                                                 std::get<2>(currentEC),
                                                 EulerCromerPend.get_nrg()});
    }

    std::ofstream Outfile(FileName);
    if (Outfile.fail())
    {
        std::cerr << "Failed to open " << FileName << ". Check that you have "
                  << "the necessary permisions and that the path exists."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    Outfile << "time,e_pos,e_nrg,ec_pos,ec_nrg"
            << std::endl;
    for (auto r = Results.begin(); r != Results.end(); r++)
    {
        for (auto a = (*r).begin(); a != (*r).end(); a++)
        {
            Outfile << *a << ",";
        }
        Outfile << std::endl;
    }
    Outfile.close();

    return 0;
}
