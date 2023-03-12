#ifndef WILL_EULERCROMERPENDULUM_H
#define WILL_EULERCROMERPENDULUM_H

#include <pendulum.h>

namespace will {

class EulerCromerPendulum : public Pendulum
{
public:
    using Pendulum::Pendulum;

    std::tuple<double, double, double> advance(double step) override;

private:
    // intentionally blank
};

} // namespace will

#endif // WILL_EULERCROMERPENDULUM_H
