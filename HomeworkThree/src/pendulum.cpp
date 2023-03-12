#include <pendulum.h>

using namespace will;

Pendulum::Pendulum(std::pair<double, double> init_cond,
                      double init_step)
{
    assert(init_step != 0.0);

    auto[init_vel, init_pos] = init_cond;
    State = {
        accel(init_pos),
        init_vel,
        init_pos
    };

    step_size = init_step;

    energy = 0.5 * std::pow(L, 2) * std::pow(init_vel, 2) -
             g * L * std::cos(init_pos);
}

std::tuple<double, double, double> Pendulum::advance_steps(int num_steps)
{
    assert(num_steps != 0);

    for (int i = 0; i < num_steps; i++)
    {
        advance(step_size);
    }

    return State;
}

std::tuple<double, double, double> Pendulum::get_state()
{
    return State;
}

double Pendulum::get_nrg()
{
    return energy;
}

double Pendulum::accel(double position)
{
    return(
        -(g / L) * std::sin(position)
    );
}

double Pendulum::vel(double vel)
{
    return vel;
}
