import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Data = pd.read_csv(
    "../data/signal.txt",
    header=None, names = ["x", "orig", "noisy"],
    dtype = {"x":np.float64, "orig":np.float64, "noisy":np.float64}
)

print(Data["x"])

fig, ax = plt.subplots()
ax.plot(Data["x"], Data["orig"])
