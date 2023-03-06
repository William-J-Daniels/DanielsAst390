#include <ode_integrator.h>

using namespace will;

OdeIntegrator::OdeIntegrator(std::function<double(double)> init_DiffEq,
                             std::pair<double, double> init_cond,
                             double init_step)
{
    assert(step_size != 0.0);

    DiffEq = init_DiffEq;
    state = init_cond;
    step_size = init_step;
}

std::pair<double, double> OdeIntegrator::advance_steps(int num_steps)
{
    assert(num_steps != 0);

    for (int i = 0; i < num_steps; i++)
    {
        advance(step_size);
    }

    return state;
}

void OdeIntegrator::set_state(std::pair<double, double> new_state)
{
    state = new_state;
}

const std::pair<double, double> OdeIntegrator::get_state()
{
    return state;
}
