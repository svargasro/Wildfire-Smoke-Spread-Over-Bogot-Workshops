import numpy as np
import os
import sys
from math import sqrt
import subprocess

#-------------------------------CONSTANTES GLOBALES------------------------
# Dimensiones de la cuadrícula
Lx = 100  # 100
Ly = int(Lx * 1.4)  # 140
iter_per_hour = 3
# iter_per_hour = 9604

Q = 9  # Número de direcciones en el espacio de velocidades

t_hour = 120
# t_hour = 240
Ux = np.zeros(Lx * Ly * t_hour)
Uy = np.zeros(Lx * Ly * t_hour)  # Velocidades en x y y

# Variables para añadir fuentes del incendio
id_list = []
rho_f = []
index_f = 0
IsData = True

C = 1.0
Cs = C / sqrt(3)  # Velocidad de la onda sonora
Cs2 = Cs ** 2  # Velocidad de la onda sonora al cuadrado

#-------------------------------VARIABLES GLOBALES------------------------
D = 0.016  # Coeficiente de difusión
tau = None  # Tiempo de relajación
Utau = None  # Inverso del tiempo de relajación
UmUtau = None  # 1 - 1/tau

# Variable para contar iteraciones por día
iteraciones_por_dia = 0

#-------------------------------CLASE LatticeBoltzmann------------------------
class LatticeBoltzmann:
    def __init__(self):
        self.w = np.zeros(Q)  # Pesos de las direcciones
        self.Vx = np.zeros(Q, dtype=int)
        self.Vy = np.zeros(Q, dtype=int)
        self.f = np.zeros(Lx * Ly * Q)
        self.fnew = np.zeros(Lx * Ly * Q)

        # Set the weights
        self.w[0] = 4.0 / 9
        self.w[1:5] = 1.0 / 9
        self.w[5:9] = 1.0 / 36

        # Set the velocity vectors
        self.Vx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
        self.Vy = [0, 0, 1, 0,-1, 1, 1, -1, -1]

    def n(self, ix, iy, i):
        return (ix * Ly + iy) * Q + i

    def rho(self, ix, iy, use_new):
        total_sum = 0.0
        for i in range(self.Q):
            n0 = self.n(ix, iy, i)
            if use_new:
                total_sum += self.fnew[n0]
            else:
                total_sum += self.f[n0]
        return total_sum

    def Jx(self, ix, iy, use_new):
        total_sum = 0.0
        for i in range(self.Q):
            n0 = self.n(ix, iy, i)
            if use_new:
                total_sum += self.Vx[i] * self.fnew[n0]
            else:
                total_sum += self.Vx[i] * self.f[n0]
        return total_sum

    def Jy(self, ix, iy, use_new):
        total_sum = 0.0
        for i in range(self.Q):
            n0 = self.n(ix, iy, i)
            if use_new:
                total_sum += self.Vy[i] * self.fnew[n0]
            else:
                total_sum += self.Vy[i] * self.f[n0]
        return total_sum

    def feq(self, rho0, Ux0, Uy0, i):
        UdotVi = Ux0 * self.Vx[i] + Uy0 * self.Vy[i]
        U2 = Ux0 * Ux0 + Uy0 * Uy0
        result = (rho0 * self.w[i] *
                  (1 + UdotVi / self.Cs2 +
                   (UdotVi * UdotVi) / (2.0 * self.Cs2 * self.Cs2) -
                   U2 / (2.0 * self.Cs2)))
        return result

    def Start(self, rho0, Ux0, Uy0):
        for ix in range(Lx):
            for iy in range(Ly):
                for i in range(self.Q):
                    n0 = self.n(ix, iy, i)
                    self.f[n0] = self.feq(rho0, Ux0, Uy0, i)

    def Collision(self):
        for ix in range(Lx):
            for iy in range(Ly):
                rho0 = self.rho(ix, iy, False)
                Ux0 = self.Jx(ix, iy, False) / rho0
                Uy0 = self.Jy(ix, iy, False) / rho0
                
                for i in range(self.Q):
                    n0 = self.n(ix, iy, i)
                    self.fnew[n0] = self.UmUtau * self.f[n0] + self.Utau * self.feq(rho0, Ux0, Uy0, i)

    def ImposeFields(self, t):
        auxT = (t+1)//self.iter_per_hour
        for ix in range(Lx):
            for iy in range(Ly):
                index = ix*Ly+iy + Lx*Ly*auxT
                Ux0 = self.Ux[index]
                Uy0 = self.Uy[index]
                rho0 = self.rho(ix, iy, True)
                
                for i in range(self.Q):
                    n0 = self.n(ix, iy, i)
                    if ((ix*Ly)+iy) == self.id_list[index]:
                        self.fnew[n0] = self.feq(self.rho_f[index], Ux0, Uy0, i)
                    else:
                        self.fnew[n0] = self.feq(rho0, Ux0, Uy0, i)

    def Advection(self):
        for ix in range(Lx):
            for iy in range(Ly):
                for i in range(self.Q):
                    ixnext = ix + self.Vx[i]
                    iynext = iy + self.Vy[i]
                    if 0 <= ixnext < Lx and 0 <= iynext < Ly:
                        n0 = self.n(ix, iy, i)
                        n0next = self.n(ixnext, iynext, i)
                        self.f[n0next] = self.fnew[n0]

    def PrintData(self, NameFile, t):
        with open(NameFile, 'w') as MyFile:
            step = 1
            for ix in range(0, Lx, step):
                for iy in range(0, Ly, step):
                    rho0 = self.rho(ix, iy, False)
                    MyFile.write(f"{ix} {iy} {rho0}\n")

    def PrintFrame(self, t):
        import os
        import subprocess
        
        if not os.path.exists('frames'):
            os.makedirs('frames')
        
        with open("frame_script.gp", "w") as GnuplotScript:
            GnuplotScript.write("""
set terminal pngcairo size 800,800 enhanced font 'Verdana,10'
set output 'frames/density_"""+str(t).zfill(3)+""".png'
set pm3d map
set size ratio -1
set xrange [0:""+str(Ly)+"]""")
            
            GnuplotScript.write(f"set yrange [0:{Lx}]")
            GnuplotScript.write("\nset cbrange [0:5000000]")
            GnuplotScript.write("\nset palette defined (0 'black', 1'red', 2 'orange', 3 'yellow', 4 'white')")
            GnuplotScript.write(f"\nset title 'Densidad en t = "+str(t)+"'")
            GnuplotScript.write('\nplot "data/density_' + str(t).zfill(3) + '.dat" u 2:1:3 w image')

        subprocess.run(["gnuplot", "frame_script.gp"])


