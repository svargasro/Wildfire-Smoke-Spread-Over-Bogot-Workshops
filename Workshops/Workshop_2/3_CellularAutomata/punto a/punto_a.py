import numpy as np
import matplotlib.pyplot as plt

Target= 'punto_a.txt'

time = np.loadtxt(Target)[:,0]#[:,0] para seleccionar la columna 1
sigma_x = np.loadtxt(Target)[:,1]
sigma_y = np.loadtxt(Target)[:,2]
sigma2 = np.loadtxt(Target)[:,3]

coeff_x = np.polyfit(time, sigma_x, 1)
line_x = np.polyval(coeff_x, time)
coeff_y = np.polyfit(time, sigma_y, 1)
line_y = np.polyval(coeff_y, time)

coeff_2 = np.polyfit(time, sigma2, 1)
line_2 = np.polyval(coeff_2, time)

plt.style.use('ggplot')

plt.title(r'Varianza en x e y vs t', fontsize=20)
plt.xlabel('Tiempo (t)')
plt.ylabel(r'Varianza ($\sigma^{2}_x$), ($\sigma^{2}_y$)')
plt.plot(time, sigma_x, color='#EA5910', label=r'$\sigma^{2}_x$')
plt.plot(time, sigma_y, color='#156C18', label=r'$\sigma^{2}_y$')
plt.plot(time, line_x, color='#993A0B', linestyle='--', label=r'Ajuste lineal para $\sigma^2_x$, m = {:.2f}'.format(coeff_x[0]))
plt.plot(time, line_y, color='#0A440C', linestyle='--', label=r'Ajuste lineal para $\sigma^2_y$, m = {:.2f}'.format(coeff_y[0]))
plt.legend(facecolor='white', edgecolor='black', fontsize='small')
plt.savefig(Target.replace('.txt', '_xy.png'))
plt.close()

plt.title(r'Varianza vs t', fontsize=20)
plt.xlabel('Tiempo (t)')
plt.ylabel(r'Varianza ($\sigma^{2}$)')
plt.plot(time, sigma2, color='#100E53', label=r'$\sigma^{2}_x$')
plt.plot(time, line_2, color='black', linestyle='--', label=r'Ajuste lineal para $\sigma^2$, m = {:.2f}'.format(coeff_2[0]))
plt.legend(facecolor='white', edgecolor='black', fontsize='small')
plt.savefig(Target.replace('.txt', '.png'))
plt.close()



