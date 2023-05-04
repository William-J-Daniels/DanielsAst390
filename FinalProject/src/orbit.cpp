#include <orbit.h>
#include <iostream>

using namespace will;

Orbit::Orbit(double a, double e) : a(a), e(e)
{
    /* Sets the semi major axis and eccentricity and uses those values to
     * compute initial conditons for the integrators.
     */

    init_pos = a * (1.0 - e);
    init_vel = -std::sqrt((G / a) * (1.0 + e) / (1.0 - e));

    S.xpos = 0.0;
    S.ypos = init_pos;

    S.xvel = init_vel;
    S.yvel = 0.0;

    S.time = 0.0;
    S.energy = energy();
}

std::vector<State> Orbit::EulerCromer(unsigned periods, unsigned divisions)
{
    auto History = std::vector<State>(periods*divisions + 2);
    History[0] = S;

    const double dt = std::pow(a, 1.5) / divisions;

    for (unsigned i = 1; i <= periods*divisions + 1; i++)
    {
        S.xvel += dt * acceleration(S.xpos);
        S.yvel += dt * acceleration(S.ypos);

        S.xpos += dt * S.xvel;
        S.ypos += dt * S.yvel;

        S.time = dt * i;
        S.energy = energy();

        History[i] = S;
    }

    return History;
}

std::vector<State> Orbit::Verlet(unsigned periods, unsigned divisions)
{
    auto History = std::vector<State>(periods*divisions + 2);
    History[0] = S;

    const double dt = std::pow(a, 1.5) / divisions;

    for (unsigned i = 1; i <= periods*divisions + 1; i++)
    {
        S.xvel += 0.5*dt * acceleration(S.xpos);
        S.yvel += 0.5*dt * acceleration(S.ypos);

        S.xpos += dt * S.xvel;
        S.ypos += dt * S.yvel;

        S.xvel += 0.5*dt * acceleration(S.xpos);
        S.yvel += 0.5*dt * acceleration(S.ypos);

        S.time = dt * i;
        S.energy = energy();

        History[i] = S;
    }

    return History;
}

std::vector<State> Orbit::Yoshida(unsigned periods, unsigned divisions)
{
    auto History = std::vector<State>(periods*divisions + 2);
    History[0] = S;

    const double dt = std::pow(a, 1.5) / divisions;

    // coefficients derived by Yoshida
    double c1 = 0.6756, c2 = -0.1756,
           d1 = 1.3512, d2 = -1.7024;

    for (unsigned i = 1; i <= periods*divisions + 1; i++)
    {
        // step 1
        S.xpos += c1 * S.xvel * dt;
        S.ypos += c1 * S.yvel * dt;

        S.xvel += d1 * acceleration(S.xpos) * dt;
        S.yvel += d1 * acceleration(S.ypos) * dt;

        // step 2
        S.xpos += c2 * S.xvel * dt;
        S.ypos += c2 * S.yvel * dt;

        S.xvel += d2 * acceleration(S.xpos) * dt;
        S.yvel += d2 * acceleration(S.ypos) * dt;

        // step 3
        S.xpos += c2 * S.xvel * dt;
        S.ypos += c2 * S.yvel * dt;

        S.xvel += d1 * acceleration(S.xpos) * dt;
        S.yvel += d1 * acceleration(S.ypos) * dt;

        // final step
        S.xpos += c1 * S.xvel * dt;
        S.ypos += c1 * S.yvel * dt;
        // no velocity update

        S.time = dt * i;
        S.energy = energy();

        History[i] = S;
    }

    return History;
}

void Orbit::reset()
{
    S.xpos = 0.0;
    S.ypos = init_pos;

    S.xvel = init_vel;
    S.yvel = 0.0;

    S.time = 0.0;
    S.energy = energy();
}

double Orbit::acceleration(double pos)
{
    return(-G * pos / std::pow(S.xpos*S.xpos + S.ypos*S.ypos, 3.0/2.0));
}

double Orbit::energy()
{
    return(
        0.5 * (S.xvel*S.xvel + S.yvel*S.yvel) -      // kinetic
        G / std::sqrt(S.xpos*S.xpos + S.ypos*S.ypos) // potential
    );
}
