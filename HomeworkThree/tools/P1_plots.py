import pandas as pd
import matplotlib.pyplot as plt


Data = pd.read_csv("../data/P1.csv", index_col=False)


fig, axs = plt.subplots(2, 1)

axs[0].set_title("Displacement and energy of Euler pendulum")
axs[0].set_ylabel("Displacement [rad]")
axs[0].set_xlabel("Time [s]")

axs[0].plot(Data["time"], Data["e_pos"],
            color="black", linestyle='-',
            label="Displacement")

ax0 = axs[0].twinx()
ax0.plot(Data["time"], Data["e_nrg"],
         color="black", linestyle=":",
         label="Energy")

lines, labels = axs[0].get_legend_handles_labels()
lines2, labels2 = ax0.get_legend_handles_labels()
axs[0].legend(lines + lines2, labels + labels2, loc="upper right",
              facecolor='white', framealpha=1)
# http://stackoverflow.com/a/10129461/1319447


axs[1].set_title("Displacement and energy of Euler-Cromer pendulum")
axs[1].set_ylabel("Displacement [rad]")
axs[1].set_xlabel("Time [s]")

axs[1].plot(Data["time"], Data["ec_pos"],
            color="black", linestyle="-",
            label="Displacement")

ax1 = axs[1].twinx()
ax1.plot(Data["time"], Data["ec_nrg"],
         color="black", linestyle = ":",
         label="Energy")

lines3, labels3 = axs[1].get_legend_handles_labels()
lines4, labels4 = ax1.get_legend_handles_labels()
axs[1].legend(lines3 + lines4, labels3 + labels4,
              loc="upper right", facecolor='white', framealpha=1)


fig.tight_layout()
fig.savefig("../data/P1_pos.jpeg")
