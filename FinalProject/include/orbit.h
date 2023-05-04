#ifndef WILL_ORBIT
#define WILL_ORBIT

#include <cmath>
#include <vector>

namespace will {

struct State
{
    /* A simple struct to keep code organized by holding values needed to
     * perform integration.
     */
    double xpos, ypos, xvel, yvel,
           time, energy;
};

class Orbit
{
    /* A class describing the orbit of a planet around a star of one solar mass.
     * Provides the following integrators:
     * - EulerCromer (1st order)
     * - Verlet      (2nd order)
     * - Yoshida     (4th order)
     * and methods to assess their stability.
     */
public:
    // constructors
    Orbit() = delete;
    Orbit(double a, double e);

    // integrators
    std::vector<State> EulerCromer(unsigned periods, unsigned divisions);
    std::vector<State> Verlet(unsigned periods, unsigned divisions);
    std::vector<State> Yoshida4(unsigned periods, unsigned divisions);

    // utility
    void reset();
    double ic_error();

private:
    static constexpr double G = 4.0 * M_PI*M_PI;
    double a, e; // semi-major axis and eccentricity
    double init_pos, init_vel; // initial conditions

    State S; // struct for state variables

    // helper functions
    double acceleration(double pos);
    double energy();
};

} // namespace will

#endif
