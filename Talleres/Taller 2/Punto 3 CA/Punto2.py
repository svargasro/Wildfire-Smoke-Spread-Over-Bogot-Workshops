import numpy as np
import matplotlib.pyplot as plt

time = np.loadtxt('Punto2.txt')[:,0]#[:,0] para seleccionar la columna 1
sigma_x = np.loadtxt('Punto2.txt')[:,1]
sigma_y = np.loadtxt('Punto2.txt')[:,2]
sigma2 = np.loadtxt('Punto2.txt')[:,3]

plt.style.use('ggplot')

plt.title(r'Varianza en x e y vs t', fontsize=20)
plt.xlabel('Tiempo (t)')
plt.ylabel(r'Varianza ($\sigma^{2}_x$), ($\sigma^{2}_y$)')
plt.plot(time,sigma_x,color='#EA5910',label=r'$\sigma^{2}_x$')
plt.plot(time,sigma_y,color='#156C18',label=r'$\sigma^{2}_y$')
plt.legend(facecolor='white', edgecolor='black')
plt.savefig('Punto2_xy.png')
plt.close()

plt.title(r'Varianza vs t', fontsize=20)
plt.xlabel('Tiempo (t)')
plt.ylabel(r'Varianza ($\sigma^{2}$)')
plt.plot(time,sigma2,color='#100E53',label=r'$\sigma^{2}_x$')
plt.savefig('Punto2.png')
plt.close()



