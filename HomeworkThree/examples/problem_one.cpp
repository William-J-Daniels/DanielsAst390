#include <eulerpendulum.h>
#include <eulercromerpendulum.h>
#include <iostream>
#include <list>
#include <array>
#include <fstream>
#include <string>


int main()
{
    const std::pair<double, double> INITCOND1 = {0.0, (10.0*M_PI/180.0)};
    const std::pair<double, double> INITCOND2 = {0.0, (100.0*M_PI/180.0)};
    const double TIME_STEP = 0.1;

    auto EPend1  = will::EulerPendulum(INITCOND1,
                                       TIME_STEP);
    auto ECPend1 = will::EulerCromerPendulum(INITCOND1,
                                             TIME_STEP);
    auto EPend2  = will::EulerPendulum(INITCOND2,
                                       TIME_STEP);
    auto ECPend2 = will::EulerCromerPendulum(INITCOND2,
                                             TIME_STEP);

    std::list<std::array<double, 5>> Results1;
    std::list<std::array<double, 5>> Results2;// linked lists to store results
    // time, euler theta, euler nrg, cromer theta, cromer nrg
    // limited by IO, so prefer constant pesh_back over constant random access
    const std::string TenDeg     = "../../HomeworkThree/data/TenDeg.csv";
    const std::string HundredDeg = "../../HomeworkThree/data/HundredDeg.csv";
    // from directory of executable

    /* vectors for the variables I names QQQ1, QQQ2 would have been better and
     * easier to generalize,but we only need them twice. Maybe it will be fixed
     * later
     */

    auto CurrentE1  = EPend1.get_state();
    auto currentEC1 = ECPend1.get_state();

    auto CurrentE2  = EPend1.get_state();
    auto currentEC2 = ECPend1.get_state();

    for (int i = 0; i < 100; i++)
    {
        CurrentE1  = EPend1.advance(TIME_STEP);
        currentEC1 = ECPend1.advance(TIME_STEP);
        Results1.push_back(std::array<double, 5> {i*TIME_STEP,
                                                 std::get<2>(CurrentE1),
                                                 EPend1.get_nrg(),
                                                 std::get<2>(currentEC1),
                                                 ECPend1.get_nrg()});

        CurrentE2  = EPend2.advance(TIME_STEP);
        currentEC2 = ECPend2.advance(TIME_STEP);
        Results2.push_back(std::array<double, 5> {i*TIME_STEP,
                                                 std::get<2>(CurrentE2),
                                                 EPend2.get_nrg(),
                                                 std::get<2>(currentEC2),
                                                 ECPend2.get_nrg()});
    }

    std::ofstream TenOut(TenDeg);
    if (TenOut.fail())
    {
        std::cerr << "Failed to open " << TenDeg << ". Check that you have "
                  << "the necessary permisions and that the path exists."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    TenOut << "time,e_pos,e_nrg,ec_pos,ec_nrg"
           << std::endl;
    for (auto r = Results1.begin(); r != Results1.end(); r++)
    {
        for (auto a = (*r).begin(); a != (*r).end(); a++)
        {
            TenOut << *a << ",";
        }
        TenOut << std::endl;
    }
    TenOut.close();

    std::ofstream HundredOut(HundredDeg);
    if (HundredOut.fail())
    {
        std::cerr << "Failed to open " << HundredDeg << ". Check that you have "
                  << "the necessary permisions and that the path exists."
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    HundredOut << "time,e_pos,e_nrg,ec_pos,ec_nrg"
               << std::endl;
    for (auto r = Results2.begin(); r != Results2.end(); r++)
    {
        for (auto a = (*r).begin(); a != (*r).end(); a++)
        {
            HundredOut << *a << ",";
        }
        HundredOut << std::endl;
    }
    HundredOut.close();

    return 0;
}
