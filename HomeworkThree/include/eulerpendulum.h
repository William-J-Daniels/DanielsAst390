#ifndef WILL_EULERPENDULUM_H
#define WILL_EULERPENDULUM_H

#include <pendulum.h>

namespace will {

class EulerPendulum : public Pendulum
{
public:
    using Pendulum::Pendulum;

    std::tuple<double, double, double> advance(double step) override;

private:
    // intentionally blank
};

} // namespace will

#endif // WILL_EULERPENDULUM_H
