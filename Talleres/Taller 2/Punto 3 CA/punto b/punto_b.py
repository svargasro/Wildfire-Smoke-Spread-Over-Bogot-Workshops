import numpy as np
import numpy as np
import matplotlib.pyplot as plt

Target= 'punto_b.txt'

time = np.loadtxt(Target)[:,0]#[:,0] para seleccionar la columna 1
sigma_x = np.loadtxt(Target)[:,1]
sigma_y = np.loadtxt(Target)[:,2]
sigma2 = np.loadtxt(Target)[:,3]

# Calculate linear regression line for sigma_x
coeff_x = np.polyfit(time, sigma_x, 1)# polyfit(x, y, grado del polinomio)
line_x = np.polyval(coeff_x, time)# polyval(p, x) p es el polinomio y x es el valor a evaluar

# Calculate linear regression line for sigma_y
coeff_y = np.polyfit(time, sigma_y, 1)
line_y = np.polyval(coeff_y, time)

# Calculate linear regression line for sigma2
coeff_2 = np.polyfit(time, sigma2, 1)
line_2 = np.polyval(coeff_2, time)

print()

plt.style.use('ggplot')

plt.title(r'Varianza en x e y vs t', fontsize=20)
plt.xlabel('Tiempo (t)')
plt.ylabel(r'Varianza ($\sigma^{2}_x$), ($\sigma^{2}_y$)')
plt.plot(time, sigma_x, color='#EA5910', label=r'$\sigma^{2}_x$')
plt.plot(time, sigma_y, color='#156C18', label=r'$\sigma^{2}_y$')
plt.plot(time, line_x, color='black', linestyle='--', label=r'Ajuste lineal $\sigma^2_x$ con pendiente = {:.2f}'.format(coeff_x[0]))
plt.plot(time, line_y, color='black', linestyle='--', label=r'Ajuste lineal $\sigma^2_y$ con pendiente = {:.2f}'.format(coeff_y[0]))
plt.legend(facecolor='white', edgecolor='black')
plt.savefig(Target.replace('.txt', '_xy.png'))
plt.close()

plt.title(r'Varianza vs t', fontsize=20)
plt.xlabel('Tiempo (t)')
plt.ylabel(r'Varianza ($\sigma^{2}$)')
plt.plot(time, sigma2, color='#100E53', label=r'$\sigma^{2}_x$')
plt.plot(time, line_2, color='black', linestyle='--', label=r'Ajuste lineal $\sigma^2$ con pendiente = {:.2f}'.format(coeff_2[0]))
plt.legend(facecolor='white', edgecolor='black')
plt.savefig(Target.replace('.txt', '.png'))
plt.close()


