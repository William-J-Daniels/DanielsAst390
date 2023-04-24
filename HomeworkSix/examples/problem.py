import numpy as np
import matplotlib.pyplot as plt

class FDGrid:
    """a finite-difference grid"""

    def __init__(self, nx, ng=1, xmin=0.0, xmax=1.0):
        """create a grid with nx points, ng ghost points (on each end)
        that runs from [xmin, xmax]"""

        self.xmin = xmin
        self.xmax = xmax
        self.ng = ng
        self.nx = nx

        # python is zero-based.  Make easy integers to know where the
        # real data lives
        self.ilo = ng
        self.ihi = ng+nx-1

        # physical coords
        self.dx = (xmax - xmin)/(nx-1)
        self.x = xmin + (np.arange(nx+2*ng)-ng)*self.dx

        # storage for the solution
        self.a = np.zeros((nx+2*ng), dtype=np.float64)
        self.b = np.zeros((nx+2*ng), dtype=np.float64)
        self.init = np.zeros((nx+2*ng), dtype=np.float64)

    def scratch_array(self):
        """ return a scratch array dimensioned for our grid """
        return np.zeros((self.nx+2*self.ng), dtype=np.float64)

    def fill_BCs(self):
        """ fill the ghostcells with outflow boundary conditions """
        self.a[0:self.ilo] = self.a[self.ilo]
        self.a[self.ihi:-1] = self.a[self.ihi]

        self.b[0:self.ilo] = self.b[self.ilo]
        self.b[self.ihi:-1] = self.b[self.ihi]

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(self.x[self.ilo:self.ihi+1], self.init[self.ilo:self.ihi+1],
                label="Initial conditions")
        ax.plot(self.x[self.ilo:self.ihi+1], self.a[self.ilo:self.ihi+1],
                label="Method a")
        ax.plot(self.x[self.ilo:self.ihi+1], self.b[self.ilo:self.ihi+1],
                label="Method b")
        ax.legend()
        return fig

def upwind_advection(nx, C_init, tmax=1.0, init_cond=None):
    """solve the linear advection equation using FTCS.  You are required
    to pass in a function f(g), where g is a FDGrid object that sets up
    the initial conditions"""

    g = FDGrid(nx)

    # initialize the data
    init_cond(g)

    g.init[:] = g.a[:]

    # evolution loop
    anew = g.scratch_array()
    bnew = g.scratch_array()

    # while loop for a
    t = 0.0
    C = C_init
    while t < tmax:
        dt = C*g.dx/g.a.max()

        if t + dt > tmax:
            dt = tmax - t
            C = g.a.max()*dt/g.dx

        # fill the boundary conditions
        g.fill_BCs()

        # loop over zones: note since we are periodic and both endpoints
        # are on the computational domain boundary, we don't have to
        # update both g.ilo and g.ihi -- we could set them equal instead.
        # But this is more general
        for i in range(g.ilo, g.ihi+1):
            anew[i] = g.a[i] - (dt/g.dx)*g.a[i]*(g.a[i] - g.a[i-1]) # mod a

        # store the updated solution
        g.a[:] = anew[:]

        t += dt

    # while loop for b
    t = 0.0
    C = C_init
    while t < tmax:
        # print("inside problem loop")
        dt = C*g.dx/g.b.max()

        if t + dt > tmax:
            dt = tmax - t
            C = g.b.max()*dt/g.dx

        # fill the boundary conditions
        g.fill_BCs()

        # loop over zones: note since we are periodic and both endpoints
        # are on the computational domain boundary, we don't have to
        # update both g.ilo and g.ihi -- we could set them equal instead.
        # But this is more general
        for i in range(g.ilo, g.ihi+1):
            bnew[i] = g.b[i] - (0.5*dt/g.dx)*(g.b[i]**2 - g.b[i-1]**2) # mod b

        # store the updated solution
        g.b[:] = bnew[:]

        t += dt

    return g

def shock(g):
    g.a[:] = 2.0
    g.a[g.x > 0.5] = 1.0

    g.b[:] = 2.0
    g.b[g.x > 0.5] = 1.0

C = 0.5
tmax = 0.1

for nx in [32, 64, 128, 256, 512]: # 32, 64, 128, 256, 512
    g = upwind_advection(nx, C, tmax=tmax, init_cond=shock)
    Sa = (g.x[g.a <= 1.5][0] - g.x[g.init <= 1.5][0]) / tmax
    Sb = (g.x[g.b <= 1.5][1] - g.x[g.init <= 1.5][0]) / tmax
    fig = g.plot()
    fig.tight_layout()
    fig.savefig("../data/"+str(nx)+".jpeg")
    print(
        "With "+str(nx)+" points, method a predicts a shock speed of "+str(Sa)+
        " and method b a shock speed of "+str(Sb)+"."
        )
