import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

class Grid:
    """ a simple class to hold cell-centered finite-difference  /
    finite-volume data """

    def __init__(self, nx, ng, xmin=0.0, xmax=1.0):

        self.xmin = xmin
        self.xmax = xmax
        self.ng = ng
        self.nx = nx

        # python is zero-based.  Make easy integers to know where the
        # real data lives
        self.ilo = ng
        self.ihi = ng+nx-1

        # physical coords -- cell-centered
        self.dx = (xmax - xmin)/(nx)
        self.x = xmin + (np.arange(nx+2*ng)-ng+0.5)*self.dx

        # storage for the solution
        self.phi = np.zeros((nx+2*ng), dtype=np.float64)

    def scratch_array(self):
        """ return a scratch array dimensioned for our grid """
        return np.zeros((self.nx+2*self.ng), dtype=np.float64)

    def fill_BCs(self):
        """ fill the ghostcells with zero gradient (Neumann)
            boundary conditions """
        self.phi[0:self.ilo] = self.phi[self.ilo]
        self.phi[self.ihi+1:] = self.phi[self.ihi]

    def norm(self, e):
        """ return the norm of quantity e which lives on the grid """
        if not len(e) == (2*self.ng + self.nx):
            return None

        return np.sqrt(self.dx*np.sum(e[self.ilo:self.ihi+1]**2))

def gaussian_ic(g, k, t=0.0, t0=1.e-4, phi1=1.0, phi2=2.0):
    xc = 0.5*(g.xmin + g.xmax)
    return (phi2 - phi1) * (
        np.sqrt(t0/(t + t0)) * np.exp(-0.25 * (g.x - xc)**2 / (k * (t + t0)))
        ) + phi1

def implicit_step(gr, phi, k, dt):
    """ diffuse phi implicitly through timestep dt """

    phinew = gr.scratch_array()

    alpha = k * dt / gr.dx**2

    # create the RHS of the matrix
    R = phi[gr.ilo:gr.ihi+1]

    # create the diagonal, d+1 and d-1 parts of the matrix
    d = (1.0 + 2.0*alpha)*np.ones(gr.nx)

    u = -alpha*np.ones(gr.nx)
    u[0] = 0.0

    l = -alpha*np.ones(gr.nx)
    l[-1] = 0.0

    # set the boundary conditions by changing the matrix elements

    # homogeneous neumann
    d[0] = 1.0 + alpha
    d[-1] = 1.0 + alpha

    # solve
    A = np.matrix([u, d, l])

    phinew[gr.ilo:gr.ihi+1] = linalg.solve_banded((1, 1), A, R)

    return phinew

def diffuse_implicit(nx, k, C, tmax, init_cond):
    """
    the main evolution loop.  Evolve

     phi_t = k phi_{xx}

    from t = 0 to tmax
    """

    # create the grid
    ng = 0

    g = Grid(nx, ng)

    # time info
    dt = C * 0.5 *g.dx**2 / k
    t = 0.0

    # initialize the data
    g.phi[:] = init_cond(g, k)

    while t < tmax:

        # make sure we end right at tmax
        if t + dt > tmax:
            dt = tmax - t

        # diffuse for dt
        phinew = implicit_step(g, g.phi, k, dt)

        g.phi[:] = phinew[:]
        t += dt

    return g

################################################################################

C = 0.8
nx = 64
k = 1
t_diffuse = (1.0/nx)**2 / k

fig, ax = plt.subplots()
tmax = 10 * t_diffuse

Cs = [512 ,256, 128, 64, 32, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]
Ns = [16, 32, 64, 128, 256]
dts = []
errors = []

for C in Cs:
# for nx in Ns:
    g = diffuse_implicit(nx, k, C, tmax, gaussian_ic)
    phi_analytic = gaussian_ic(g, k, t=tmax)
    dts.append(C * 0.5 *g.dx**2 / k)
    # dts.append(g.dx)
    errors.append(g.norm(g.phi - phi_analytic))

ax.scatter(dts, errors, marker="x", label="original")

# redifinition of implicit step function
def implicit_step(gr, phi, k, dt):
    """ diffuse phi implicitly through timestep dt """

    phinew = gr.scratch_array()

    alpha = k * dt / gr.dx**2

    # create the RHS of the matrix
    R = phi[gr.ilo:gr.ihi+1]

    # create the diagonal, d+1 and d-1 parts of the matrix
    d = (1.0 + alpha)*np.ones(gr.nx)

    u = -0.5*alpha*np.ones(gr.nx)
    u[0] = 0.0

    l = -0.5*alpha*np.ones(gr.nx)
    l[-1] = 0.0

    # set the boundary conditions by changing the matrix elements

    # homogeneous neumann
    d[0] = 1.0 + 0.5*alpha
    d[-1] = 1.0 + 0.5*alpha

    # solve
    A = np.matrix([u, d, l])

    phinew[gr.ilo:gr.ihi+1] = linalg.solve_banded((1, 1), A, R)

    return phinew

dts = []
errors = []

for C in Cs:
# for nx in Ns:
    g = diffuse_implicit(nx, k, C, tmax, gaussian_ic)
    phi_analytic = gaussian_ic(g, k, t=tmax)
    dts.append(C * 0.5 *g.dx**2 / k)
    # dts.append(g.dx)
    errors.append(g.norm(g.phi - phi_analytic))

ax.scatter(dts, errors, marker="x", label="new")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$\Delta t$")
ax.set_ylabel("Absolute error")

ax.legend()

fig.savefig("../data/TEST.jpeg")
