#ifndef WILL_PENDULUM_H
#define WILL_PENDULUM_H

#include <tuple>
#include <cassert>
#include <cmath>

#include <iostream>

namespace will{

class Pendulum
{
    // Absgract pendulum class
public:
    Pendulum(){ }
    Pendulum(std::pair<double, double> init_cond,
             double init_step);
    virtual ~Pendulum(){ }

    virtual std::tuple<double, double, double> advance(double step) = 0;
    std::tuple<double, double, double> advance_steps(int num_steps);

    std::tuple<double, double, double> get_state();
    double get_nrg();

protected:
    static constexpr double g = 10.0; // m s^-2
    double L = 10.0; // meters

    double step_size; // seconds

    std::tuple<double, double, double> State;
            // accel , vel   , pos
    double energy = 0.0;

    double accel(double position);
    double vel(double velocity);
};

} // namespace will

#endif // WILL_PENDULUM_H
