import numpy as np
import matplotlib.pyplot as plt


x_s, y_s, x_j, y_j, x_t, y_t, t = np.genfromtxt('data.txt', unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6))

plt.style.use('seaborn-v0_8')

fig, axes = plt.subplots(figsize=(6, 6))


axes.plot(x_t, y_t, ".", color="cyan", label=r"$(x',y')_{T}$")
axes.set_xlabel("x'", fontsize=12)
axes.set_ylabel("y'", fontsize=12)
axes.legend(loc='upper right')
axes.grid(True, linestyle='--')
axes.set_title("y' vs. x'. Perturbado en 5 partes por mil.", fontsize=14)
axes.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()


fig.savefig('plot.pdf')
