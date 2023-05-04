#include <orbit.h>
#include <iostream>
#include <fstream>
#include <limits>

const unsigned PERIODS = std::pow(2, 6), DIVISIONS = 100;

int main()
{
    auto O = will::Orbit(1.0, 0.7);

    auto EulerCromerResult = O.EulerCromer(PERIODS, DIVISIONS);

    std::ofstream EulerCromerFile("../../FinalProject/data/EulerCromer.csv");
    EulerCromerFile << "xpos,ypos,xvel,yvel,time,energy" << std::endl;
    for (auto& S : EulerCromerResult)
        EulerCromerFile << S.xpos << "," << S.ypos << "," <<
                           S.xvel << "," << S.yvel << "," <<
                           S.time << "," << S.energy << std::endl;
    EulerCromerFile.close();


    O.reset();
    auto VerletResult = O.Verlet(PERIODS, DIVISIONS);

    std::ofstream VerletFile("../../FinalProject/data/Verlet.csv");
    VerletFile << "xpos,ypos,xvel,yvel,time,energy" << std::endl;
    for (auto& S : VerletResult)
        VerletFile << S.xpos << "," << S.ypos << "," <<
                      S.xvel << "," << S.yvel << "," <<
                      S.time << "," << S.energy << std::endl;
    VerletFile.close();


    O.reset();
    auto Yoshida4Result = O.Yoshida4(PERIODS, DIVISIONS);

    std::ofstream Yoshida4File("../../FinalProject/data/Yoshida4.csv");
    Yoshida4File << "xpos,ypos,xvel,yvel,time,energy" << std::endl;
    for (auto& S : Yoshida4Result)
        Yoshida4File << S.xpos << "," << S.ypos << "," <<
                       S.xvel << "," << S.yvel << "," <<
                       S.time << "," << S.energy << std::endl;
    Yoshida4File.close();

    return 0;
}
