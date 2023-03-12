#include <velocityverletpendulum.h>

using namespace will;

std::tuple<double, double, double> VelocityVerletPendulum::advance(double step)
{
    auto[acceleration, velocity, position] = State;

    velocity += accel(position)*(step*0.5);
    position += velocity*step;
    velocity += accel(position)*(step*0.5);
    acceleration = accel(position);

    State = {acceleration, velocity, position};

    energy = 0.5 * std::pow(L, 2) * std::pow(velocity, 2) -
             g * L * std::cos(position);

    return State;
}
