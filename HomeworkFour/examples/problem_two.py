# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# function definitions
def kernal(x, sigma):
    return(
        1.0 / (sigma * np.sqrt(2.0 * np.pi)) *
        np.exp(-0.5 * (x/sigma)**2)
    )


# read and plot data
Data = pd.read_csv(
    "../data/signal.txt",
    header=None, names = ["x", "orig", "noisy"], delim_whitespace=True,
    dtype = {"x":np.float64, "orig":np.float64, "noisy":np.float64}
)

fig1, axs1 = plt.subplots(
    2, 1, layout="tight",
    sharex=True, sharey=True
)

axs1[0].set_title("Original")
axs1[0].plot(Data["x"], Data["orig"],
             color="black")

axs1[1].set_title("Noisy")
axs1[1].plot(Data["x"], Data["noisy"],
             color="black")

fig1.savefig("../data/plotted_file.jpeg")


# plot the noisy function and the kernal together
Data["kernal"] = kernal(Data["x"], 1.0)
Data["kernal"] = Data["kernal"] + Data["kernal"].values[::-1]
Data["kernal"] = Data["kernal"] / Data["kernal"].sum()

fig2, ax2 = plt.subplots()
ax2.set_title("Noisy function and kernal")
ax2.plot(Data["x"], Data["noisy"],
         color="black")
ax2.twinx().plot(Data["x"], Data["kernal"],
         color="black", ls=":")

fig2.savefig("../data/NoiseAndKernal.jpeg")


# plot their fourier transforms
Data["fft_noisy"] = np.fft.fft(Data["noisy"])
Data["fft_kernal"] = np.fft.fft(Data["kernal"])

fig3, ax3 = plt.subplots()
ax3.set_title("Noisy function and kernal FFT")
ax3.plot(Data["x"], Data["fft_noisy"],
         color="black")
ax3.twinx().plot(Data["x"], Data["fft_kernal"],
                 color="black", ls=":")

fig3.savefig("../data/NoiseAndKernalFFT.jpeg")


# convolution and denoise
Data["conv"] = Data["fft_noisy"] * Data["fft_kernal"]
Data["denoise"] = np.fft.ifft(Data["conv"])

fig4, axs4 = plt.subplots(
    2, 1, layout="tight",
    sharex=True, sharey=True
)

axs4[0].set_title("Denoised")
axs4[0].plot(Data["x"], Data["denoise"],
             color="black")

axs4[1].set_title("Original")
axs4[1].plot(Data["x"], Data["orig"],
             color="black")

fig4.savefig("../data/denoise.jpeg")
