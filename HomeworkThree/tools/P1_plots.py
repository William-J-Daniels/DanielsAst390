import pandas as pd
import matplotlib.pyplot as plt


DataTen = pd.read_csv("../data/TenDeg.csv", index_col=False)
DataHundred = pd.read_csv("../data/HundredDeg.csv", index_col=False)


fig, axs = plt.subplots(2, 2, figsize=(7.5, 5),
                        layout="tight")
# should have used loops- didn't initially intend for subplots

axs[0,0].set_title("Euler 10 deg")
axs[0,0].set_ylabel("Displacement [rad]")
axs[0,0].set_xlabel("Time [s]")

axs[0,0].plot(DataTen["time"], DataTen["e_pos"],
            color="black", linestyle='-',
            label="Displacement")

ax00 = axs[0,0].twinx()
ax00.set_ylabel("Energy [arb]")
ax00.plot(DataTen["time"], DataTen["e_nrg"],
         color="black", linestyle=":",
         label="Energy")

lines, labels = axs[0,0].get_legend_handles_labels()
lines2, labels2 = ax00.get_legend_handles_labels()
ax00.legend(lines + lines2, labels + labels2,
              loc="upper right", framealpha=1.0, facecolor="white")
# http://stackoverflow.com/a/10129461/1319447


axs[1,0].set_title("Euler-Cromer 10 deg")
axs[1,0].set_ylabel("Displacement [rad]")
axs[1,0].set_xlabel("Time [s]")

axs[1,0].plot(DataTen["time"], DataTen["ec_pos"],
            color="black", linestyle="-",
            label="Displacement")

ax10 = axs[1,0].twinx()
ax10.set_ylabel("Energy [arb]")
ax10.plot(DataTen["time"], DataTen["ec_nrg"],
         color="black", linestyle = ":",
         label="Energy")

lines3, labels3 = axs[1,0].get_legend_handles_labels()
lines4, labels4 = ax10.get_legend_handles_labels()
ax10.legend(lines3 + lines4, labels3 + labels4,
              loc="upper right", facecolor="white", framealpha=1.0)


axs[0,1].set_title("Euler 100 deg")
axs[0,1].set_ylabel("Displacement [rad]")
axs[0,1].set_xlabel("Time [s]")

axs[0,1].plot(DataHundred["time"], DataHundred["e_pos"],
            color="black", linestyle='-',
            label="Displacement")

ax01 = axs[0,1].twinx()
ax01.set_ylabel("Energy [arb]")
ax01.plot(DataHundred["time"], DataHundred["e_nrg"],
         color="black", linestyle=":",
         label="Energy")

lines5, labels5 = axs[0,1].get_legend_handles_labels()
lines6, labels6 = ax01.get_legend_handles_labels()
ax01.legend(lines5 + lines6, labels5 + labels6,
              loc="upper right", framealpha=1.0, facecolor="white")


axs[1,1].set_title("Euler-Cromer 100 deg")
axs[1,1].set_ylabel("Displacement [rad]")
axs[1,1].set_xlabel("Time [s]")

axs[1,1].plot(DataHundred["time"], DataHundred["ec_pos"],
            color="black", linestyle="-",
            label="Displacement")

ax11 = axs[1,1].twinx()
ax11.set_ylabel("Energy [arb]")
ax11.plot(DataHundred["time"], DataHundred["ec_nrg"],
         color="black", linestyle = ":",
         label="Energy")

lines7, labels7 = axs[1,1].get_legend_handles_labels()
lines8, labels8 = ax11.get_legend_handles_labels()
ax11.legend(lines7 + lines8, labels7 + labels8,
              loc="upper right", facecolor="white", framealpha=1.0)


fig.savefig("../data/P1Plots.jpeg")
