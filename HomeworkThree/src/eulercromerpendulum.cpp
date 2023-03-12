#include <eulercromerpendulum.h>

using namespace will;

std::tuple<double, double, double> EulerCromerPendulum::advance(double step)
{
    auto[acceleration, velocity, position] = State;

    velocity += accel(position)*step;
    position += velocity*step;
    acceleration = accel(position);

    State = {acceleration, velocity, position};

    energy = 0.5 * std::pow(L, 2) * std::pow(velocity, 2) -
             g * L * std::cos(position);

    return State;
}
