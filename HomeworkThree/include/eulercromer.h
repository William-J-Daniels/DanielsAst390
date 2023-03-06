#ifndef WILL_EULERCROMER_H
#define WILL_EULERCROMER_H

#include <ode_integrator.h>

namespace will {

class EulerCromer : public OdeIntegrator
{
    // A fixed-time step implimentation of Euler's method for numerical
    // solutions to ODEs
public:
    using OdeIntegrator::OdeIntegrator;

    std::pair<double, double> advance(double step) override;

private:
    // intentionally blank
};

}

#endif // WILL_EULERCROMER_H
