#include <velocityverletpendulum.h>

int main()
{
    const std::pair<double, double> INITCOND = {0.0, (33.3*M_PI/180.0)};
    double time_step = 0.5,
           init_nrg, new_nrg;

    auto Pend = will::VelocityVerletPendulum();

    for (int i = 0; i < 5; i++)
    {
        Pend = will::VelocityVerletPendulum(INITCOND,
                                            time_step/std::pow(2.0, i));
        init_nrg = Pend.get_nrg();
        Pend.advance_steps(50);
        new_nrg = Pend.get_nrg();

        std::cout << new_nrg/init_nrg << std::endl;
    }

    return 0;
}
