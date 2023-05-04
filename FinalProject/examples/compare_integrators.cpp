#include <orbit.h>
#include <iostream>
#include <fstream>

const unsigned PERIODS = 20, DIVISIONS = 100;

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
    auto YoshidaResult = O.Yoshida(PERIODS, DIVISIONS);

    std::ofstream YoshidaFile("../../FinalProject/data/Yoshida.csv");
    YoshidaFile << "xpos,ypos,xvel,yvel,time,energy" << std::endl;
    for (auto& S : YoshidaResult)
        YoshidaFile << S.xpos << "," << S.ypos << "," <<
                       S.xvel << "," << S.yvel << "," <<
                       S.time << "," << S.energy << std::endl;
    YoshidaFile.close();

    return 0;
}
