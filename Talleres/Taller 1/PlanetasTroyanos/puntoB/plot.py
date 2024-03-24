import numpy as np
import matplotlib.pyplot as plt


x_s, y_s, x_j, y_j, t = np.genfromtxt('data.txt', unpack=True, usecols=(0, 1, 2, 3, 4)) #Se extraen los datos.

plt.style.use('seaborn-v0_8')

fig, axes = plt.subplots(5,1,figsize=(7, 7))




axes[0].plot(t, x_s, '.', color='yellow', label=r"$x^{'}(t)_{S}$")
axes[1].plot(t, y_s, '.', color='yellow', label=r"$y^{'}(t)_{S}$")
axes[2].plot(t, x_j, '.', color='black', label=r"$x^{'}(t)_{J}$")
axes[3].plot(t, y_j, '.', color='black', label=r"$y^{'}(t)_{J}$")

axes[4].plot(x_s, y_s, ".", color="yellow", label=r"$(x',y')_{S}$")
axes[4].plot(x_j, y_j, ".", color="black", label=r"$(x',y')_{J}$")





# Se ajustan demás detalles del gráfico.

for i in [0,2]:
    axes[i].set_xlabel('t', fontsize=12)
    axes[i].set_ylabel('x\'(t)', fontsize=12)
    axes[i].legend(loc='upper right')
    axes[i].grid(True, linestyle='--')
    axes[i].set_title("x' vs. t'", fontsize=14)
    axes[i].ticklabel_format(style='sci', axis='x', scilimits=(0,0))

for i in [1,3]:
    axes[i].set_xlabel('t', fontsize=12)
    axes[i].set_ylabel('y\'(t)', fontsize=12)
    axes[i].legend()
    axes[i].grid(True, linestyle='--')
    axes[i].set_title("y' vs. t", fontsize=14)
    axes[i].ticklabel_format(style='sci', axis='x', scilimits=(0,0))

axes[4].set_xlabel("x'", fontsize=12)
axes[4].set_ylabel("y'", fontsize=12)
axes[4].legend(loc='center')
axes[4].grid(True, linestyle='--')
axes[4].set_title("y' vs. x'", fontsize=14)
axes[4].ticklabel_format(style='sci', axis='x', scilimits=(0,0))



plt.tight_layout()
fig.savefig('plot.pdf')
