#ifndef WILL_ODEINTEGRATOR_H
#define WILL_ODEINTEGRATOR_H

#include <functional>
#include <utility>
#include <cassert>

namespace will {

class OdeIntegrator
{
    // Abstract class for the integration of one dimensional ODEs
public:
    OdeIntegrator(std::function<double(double)> init_DiffEq,
                  std::pair<double, double> init_cond,
                  double init_step);
    virtual ~OdeIntegrator(){ }

    virtual std::pair<double, double> advance(double step) = 0;
    std::pair<double, double> advance_steps(int num_steps);

    void set_state(std::pair<double, double> new_state);
    const std::pair<double, double> get_state();

protected:
    std::function<double(double)> DiffEq;
    std::pair<double, double> state;
    double step_size; // not made const so if i have to do variable step size i
                      // can inherit from this class
};

} // namespace will

#endif // WILL_INTEGRATOR_H
