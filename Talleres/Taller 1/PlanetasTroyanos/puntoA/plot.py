import numpy as np
import matplotlib.pyplot as plt


x_s, y_s, x_j, y_j = np.genfromtxt('data.txt', unpack=True, usecols=(0, 1, 2, 3))

plt.style.use('seaborn-v0_8')

fig, axes = plt.subplots(figsize=(6, 6))

axes.plot(x_s, y_s, '.', color='yellow', label=r'$(x,y)_{S}$')
axes.plot(x_j, y_j, '.', color='black', label=r'$(x,y)_{J}$')




# Se ajustan demás detalles del gráfico.

axes.set_xlabel('x', fontsize=12)
axes.set_ylabel('y', fontsize=12)
axes.legend(loc='upper left')
axes.grid(True, linestyle='--')
axes.set_title("y vs x", fontsize=14)
plt.tight_layout()
fig.savefig('plot.png')
