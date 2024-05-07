import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import multivariate_normal


Target= 'punto_b.txt'

time = np.loadtxt(Target)[:,0]#[:,0] para seleccionar la columna 1
sigma_x = np.loadtxt(Target)[:,1]
sigma_y = np.loadtxt(Target)[:,2]
sigma2 = np.loadtxt(Target)[:,3]

# Calculo regresión lineal para sigma2
coeff_2 = np.polyfit(time, sigma2, 1)
line_2 = np.polyval(coeff_2, time)

plt.style.use('ggplot')

plt.title(r'Varianza vs t', fontsize=20)
plt.xlabel('Tiempo (t)')
plt.ylabel(r'Varianza ($\sigma^{2}$)')
plt.plot(time, sigma2, color='#7feb5e', label=r'$\sigma^{2}_x$')
plt.plot(time, line_2, color='black', linestyle='--', label=r'Ajuste lineal $\sigma^2$ con pendiente = {:.2f}'.format(coeff_2[0]))
plt.legend(facecolor='white', edgecolor='black')
plt.savefig(Target.replace('.txt', '_varianza.png'))
plt.close()

#Grafica distribuciones iniciales y finales------------------------------------------------------------------------
parameters = np.loadtxt("parametros.txt")

mu_x_i = parameters[0][0]
mu_y_i = parameters[0][1]
sigma_x_i = parameters[0][2]
sigma_y_i = parameters[0][3]

mu_x_f = parameters[1][0]
mu_y_f = parameters[1][1]
sigma_x_f = parameters[1][2]
sigma_y_f = parameters[1][3]

mean_i = [mu_x_i, mu_y_i]
cov_i = [[sigma_x_i**2, 0], [0, sigma_y_i**2]]

mean_f = [mu_x_f, mu_y_f]
cov_f = [[sigma_x_f**2, 0], [0, sigma_y_f**2]]

x = np.linspace(mu_x_f - sigma_x_f*2.5, mu_x_f + sigma_x_f*2.5, 500)
y = np.linspace(mu_y_f - sigma_y_f*2.5, mu_y_f + sigma_y_f*2.5, 500)
X, Y = np.meshgrid(x, y)

# Calcular la densidad de probabilidad en cada punto
Z_i = multivariate_normal.pdf(np.dstack([X, Y]), mean=mean_i, cov=cov_i).reshape(X.shape)
Z_f = multivariate_normal.pdf(np.dstack([X, Y]), mean=mean_f, cov=cov_f).reshape(X.shape)

fig = plt.figure(figsize=(18, 6))
# Grafico para la distribución inicial
ax1 = fig.add_subplot(1, 2, 1, projection='3d') 
surf1 = ax1.plot_surface(X, Y, Z_i, cmap=cm.RdYlGn)
ax1.set_title('Distribución Inicial')  
ax1.set_xlabel('Eje X') 
ax1.set_ylabel('Eje Y') 
fig.colorbar(surf1, shrink=0.8, aspect=8, location='left', ax=ax1)

# Grafico para la distribución final
ax2 = fig.add_subplot(1, 2, 2, projection='3d') 
surf2 = ax2.plot_surface(X, Y, Z_f, cmap=cm.RdYlGn)
ax2.set_title('Distribución Final')  
ax2.set_xlabel('Eje X')  
ax2.set_ylabel('Eje Y') 
fig.colorbar(surf2, shrink=0.8, aspect=8, location='left',ax=ax2)
fig.suptitle('Comparación de distribuciones de probabilidad', fontsize=20)
plt.savefig(Target.replace('.txt', '_distribucion.png'))
plt.close()

# Graficas 2D vista superior  --------------------------------------------------
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(1, 2, figsize=(16, 6))

# Subtrama para Z_i
cmap = 'RdYlGn'
contourf = ax[0].contourf(X, Y, Z_i, levels=50, cmap=cmap)
divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right", size="5%", pad=0.1)
fig.colorbar(contourf, cax=cax)  # Añade la colorbar
ax[0].set_title('Distribución Inicial')
ax[0].set_xlabel('Eje X')
ax[0].set_ylabel('Eje Y')

# Subtrama para Z_f
contourf = ax[1].contourf(X, Y, Z_f, levels=50, cmap=cmap)
divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right", size="5%", pad=0.1)
fig.colorbar(contourf, cax=cax)  # Añade la colorbar
ax[1].set_title('Distribución Final')
ax[1].set_xlabel('Eje X')
ax[1].set_ylabel('Eje Y')

plt.subplots_adjust(wspace=0.5)
fig.suptitle('Distribución Gaussiana Bidimensional (Vista superior)', fontsize=20)
plt.savefig(Target.replace('.txt', '_superior2D.png'))
plt.close()

# fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# # Subtrama para Z_i
# contour0 = ax[0].contourf(X, Y, Z_i, levels=50, cmap='RdYlGn')
# ax[0].contourf(X, Y, Z_i, levels=50, cmap='RdYlGn')
# ax[0].set_title('Distribución Inicial')
# ax[0].set_xlabel('Eje X')
# ax[0].set_ylabel('Eje Y')
# cbar_ax0 = fig.add_axes([0.5, 0.15, 0.03, 0.7])
# fig.colorbar(contour0, cax=cbar_ax0, orientation='vertical')

# # Subtrama para Z_f
# contour1 = ax[1].contourf(X, Y, Z_f, levels=50, cmap='RdYlGn')
# ax[1].contourf(X, Y, Z_f, levels=50, cmap='RdYlGn')
# ax[1].set_title('Distribución Final')
# ax[1].set_xlabel('Eje X')
# ax[1].set_ylabel('Eje Y')
# cbar_ax1 = fig.add_axes([0.92, 0.15, 0.03, 0.7])  
# fig.colorbar(contour1, cax=cbar_ax1, orientation='vertical')  

# plt.subplots_adjust(left=0.1, right=0.85) 
# fig.suptitle('Distribución Gaussiana Bidimensional (Vista superior)', fontsize=20)
# plt.savefig(Target.replace('.txt', '_superior2D.png'))
# plt.close()

