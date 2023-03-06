#ifndef WILL_EULERMETHOD_H
#define WILL_EULERMETHOD_H

#include <ode_integrator.h>

namespace will {

class EulerMethod : public OdeIntegrator
{
    // A fixed-time step implimentation of Euler's method for numerical
    // solutions to ODEs
public:
    using OdeIntegrator::OdeIntegrator;

    std::pair<double, double> advance(double step) override;

private:
    // intentionally left blank
};

}

#endif // WILL_EULERMETHOD_H
