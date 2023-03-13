#include <eulerpendulum.h>

using namespace will;

std::tuple<double, double, double> EulerPendulum::advance(double step)
{
    auto[acceleration, velocity, position] = State;

    position += velocity*step;
    velocity += accel(position)*step;
    acceleration = accel(position);

    State = {acceleration, velocity, position};

    energy = 0.5 * std::pow(L, 2) * std::pow(velocity, 2) +
             g * L * (1.0 - std::cos(position));

    return State;
}
