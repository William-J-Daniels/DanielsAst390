#ifndef WILL_VELOCITYVERLETPENDULUM_H
#define WILL_VELOCITYVERLETPENDULUM_H

#include <pendulum.h>

namespace will {

class VelocityVerletPendulum : public Pendulum
{
public:
    using Pendulum::Pendulum;

    std::tuple<double, double, double> advance(double step) override;

private:
    // intentionally blank
};

}

#endif // WILL_VELOCITYVERLETPENDULUM_H
