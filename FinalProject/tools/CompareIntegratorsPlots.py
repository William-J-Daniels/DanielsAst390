import pandas as pd
import matplotlib.pyplot as plt

EulerCromerData = pd.read_csv("../data/EulerCromer.csv")
VerletData = pd.read_csv("../data/Verlet.csv")
Yoshida4Data = pd.read_csv("../data/Yoshida4.csv")

fig, ax = plt.subplots()
ax.set_title("Orbits")
ax.set_xlabel("X position")
ax.set_ylabel("Y position")

ax.scatter([0], [0], marker="*", color="black", label="Star")
ax.plot(EulerCromerData["xpos"], EulerCromerData["ypos"],
        lw=0.5, label="EulerCromer")
ax.plot(VerletData["xpos"], VerletData["ypos"],
        lw=0.5, label="Verlet")
ax.plot(Yoshida4Data["xpos"], Yoshida4Data["ypos"],
        lw=0.5, label="Yoshida4")

ax.legend()
ax.set_aspect(1.0)
fig.tight_layout()

fig.savefig("../data/Orbits.jpeg")


fig, ax = plt.subplots()
ax.set_title("Energies")
ax.set_xlabel("Time")
ax.set_ylabel("Energy")

ax.plot(EulerCromerData["time"], EulerCromerData["energy"],
        lw=0.5, label="EulerCromer")
ax.plot(VerletData["time"], VerletData["energy"],
        lw=0.5, label="Verlet")
ax.plot(Yoshida4Data["time"], Yoshida4Data["energy"],
        lw=0.5, label="Yoshida4")

ax.legend()
fig.tight_layout()

fig.savefig("../data/Energies.jpeg")
