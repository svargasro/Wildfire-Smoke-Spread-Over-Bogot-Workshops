import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


x_s, y_s, x_j, y_j, x_t, y_t, t = np.genfromtxt('data.txt', unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6))

plt.style.use('seaborn-v0_8')

fig, axes = plt.subplots(figsize=(6, 6))

axes.plot(t, x_t, ".", color="cyan", label=r"$(x',y')_{T}$")
axes.set_xlabel("t", fontsize=12)
axes.set_ylabel("x'", fontsize=12)
axes.legend(loc='upper right')
axes.grid(True, linestyle='--')
axes.set_title("t vs. x'. Perturbado en 5 partes por mil.", fontsize=14)
axes.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()


peaks, _ = find_peaks(x_t) #Se encuentran los máximos locales.

tpeak = t[peaks] #Se hallan los tiempos para dichos picos.
tpeaksD = np.roll(t[peaks], 1) #Se corren los valores del arreglo en una posición.
periodos = tpeak - tpeaksD #Se hallan el tiempo que hay entre cada pico.
periodos = periodos[1:] #Se descarta el primer valor puesto que no aporta información.
periodo2 = np.mean(periodos)
print("Periodo 2: ", periodo2)


midT = int(len(t)/2)

#Se hallan los tiempos para los cuales corresponden los dos mayores maximos partiendo el intervalo de tiempos en 2.
tMax1 = t[np.argmax(x_t[:midT])]
tMax2 = t[midT-1+np.argmax(x_t[midT:])]

#Graficamente, se ve que para esta perturbación se tienen dos picos en el segundo máximo pero
#parece que la periodicidad se alcanza en el medio de dichos dos máximos. Así, se calcula el tMin entre los dos máximos.
#Como ya conocemos todos los máximos locales, basta hacer el promedio con el máximo anterior a tMax.






tolerancia = 10
indicetMax2 = np.where(np.abs(tpeak - tMax2) <= tolerancia)[0]   #Se halla el índice al cual corresponde tMax2 de todos los picos.
tMaxAnt2 = tpeak[indicetMax2[0] -1] #Se toma el valor de tiempo del máximo anterior a tMax2



tMax2 = (tMax2 + tMaxAnt2)/2.0 #Se sobreescribe tMax2.

periodo1 = tMax2 - tMax1

print("Periodo 1: ", periodo1)


fig.savefig('plot.pdf')
