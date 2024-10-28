import numpy as np
cimport numpy as np
from libc.math cimport exp, sqrt, pow, M_PI
from libc.stdio cimport FILE, fopen, fclose, fprintf

#-------------------------------CONSTANTES GLOBALES------------------------
# Dimensiones de la cuadrícula
cdef int Lx = 100
cdef int Ly = int(Lx * 1.4)
cdef int iter_per_hour = 9604

cdef int t_hour = 240
cdef np.ndarray[np.double_t, ndim=1] Ux = np.zeros(Lx * Ly * t_hour, dtype=np.double)
cdef np.ndarray[np.double_t, ndim=1] Uy = np.zeros(Lx * Ly * t_hour, dtype=np.double)

# Número de direcciones en el espacio de velocidades
cdef int Q = 9

# Velocidad del lattice y velocidad del sonido
cdef double C = 1.0            # Velocidad característica del lattice
cdef double Cs = C / sqrt(3) # Velocidad de la onda sonora
cdef double Cs2 = Cs * Cs    # Velocidad de la onda sonora al cuadrado

#-------------------------------VARIABLES GLOBALES------------------------

# Parámetros de difusión y relajación
cdef double D = 0.016  # Coeficiente de difusión
cdef double tau       # Tiempo de relajación
cdef double Utau      # Inverso del tiempo de relajación
cdef double UmUtau    # 1 - 1/tau

#-------------------------------CLASES--------------------------------------
# Clase para el método de Lattice Boltzmann
cdef class LatticeBoltzman:
    cdef:
        double w[Q]  # Pesos de las direcciones
        int Vx[Q], Vy[Q]  # Vectores de velocidad en las direcciones x e y
        double* f
        double* fnew

    def _cinit_(self):
        # Asignación de los pesos para cada dirección
        self.w[0] = 4.0 / 9
        for i in range(1, 5):
            self.w[i] = 1.0 / 9
        for i in range(5, 9):
            self.w[i] = 1.0 / 36

        # Asignación de los vectores de velocidad
        self.Vx[8] = 1
        self.Vx[1] = 1
        self.Vx[5] = 1
        self.Vx[4] = 0
        self.Vx[0] = 0
        self.Vx[2] = 0
        self.Vx[7] = -1
        self.Vx[3] = -1
        self.Vx[6] = -1

        self.Vy[8] = -1
        self.Vy[1] = 0
        self.Vy[5] = 1
        self.Vy[4] = -1
        self.Vy[0] = 0
        self.Vy[2] = 1
        self.Vy[7] = -1
        self.Vy[3] = 0
        self.Vy[6] = 1

        # Creación de los arreglos dinámicos para las funciones de distribución
        cdef int ArraySize = Lx * Ly * Q
        self.f = <double *> malloc(ArraySize * sizeof(double))
        self.fnew = <double *> malloc(ArraySize * sizeof(double))
        for i in range(ArraySize):
            self.f[i] = 0.0
            self.fnew[i] = 0.0

    def _dealloc_(self):
        free(self.f)
        free(self.fnew)

    cdef int n(self, int ix, int iy, int i):
        return (ix * Ly + iy) * Q + i

    cdef double rho(self, int ix, int iy, bint UseNew):
        cdef double sum = 0.0
        cdef int i, n0
        for i in range(Q):
            n0 = self.n(ix, iy, i)
            if UseNew:
                sum += self.fnew[n0]
            else:
                sum += self.f[n0]
        return sum

    cdef double Jx(self, int ix, int iy, bint UseNew):
        cdef double sum = 0.0
        cdef int i, n0
        for i in range(Q):
            n0 = self.n(ix, iy, i)
            if UseNew:
                sum += self.Vx[i] * self.fnew[n0]
            else:
                sum += self.Vx[i] * self.f[n0]
        return sum

    cdef double Jy(self, int ix, int iy, bint UseNew):
        cdef double sum = 0.0
        cdef int i, n0
        for i in range(Q):
            n0 = self.n(ix, iy, i)
            if UseNew:
                sum += self.Vy[i] * self.fnew[n0]
            else:
                sum += self.Vy[i] * self.f[n0]
        return sum

    cdef double feq(self, double rho0, double Ux0, double Uy0, int i):
        cdef double UdotVi = Ux0 * self.Vx[i] + Uy0 * self.Vy[i]
        cdef double U2 = Ux0 * Ux0 + Uy0 * Uy0
        return rho0 * self.w[i] * (1 + UdotVi / Cs2 + (UdotVi * UdotVi) / (2.0 * Cs2 * Cs2) - U2 / (2.0 * Cs2))

    cpdef Start(self, double rho0, double Ux0, double Uy0, double mu_x, double mu_y, double sigma_x, double sigma_y):
        cdef int ix, iy, i, n0
        cdef double gauss_x, gauss_y, rho
        for ix in range(Lx):
            for iy in range(Ly):
                gauss_x = exp(-0.5 * pow((ix - mu_x) / sigma_x, 2)) / (sigma_x * sqrt(2 * M_PI))
                gauss_y = exp(-0.5 * pow((iy - mu_y) / sigma_y, 2)) / (sigma_y * sqrt(2 * M_PI))
                rho = rho0 * gauss_x * gauss_y
                for i in range(Q):
                    n0 = self.n(ix, iy, i)
                    self.f[n0] = self.feq(rho, Ux0, Uy0, i)

    cpdef Collision(self):
        cdef int ix, iy, i, n0
        cdef double rho0, Ux0, Uy0
        for ix in range(Lx):
            for iy in range(Ly):
                rho0 = self.rho(ix, iy, False)
                Ux0 = self.Jx(ix, iy, False) / rho0
                Uy0 = self.Jy(ix, iy, False) / rho0
                for i in range(Q):
                    n0 = self.n(ix, iy, i)
                    self.fnew[n0] = UmUtau * self.f[n0] + Utau * self.feq(rho0, Ux0, Uy0, i)

    cpdef ImposeFields(self, int t):
        cdef int ix, iy, i, index, auxT = (t + 1) / iter_per_hour
        cdef double rho0, Ux0, Uy0
        for ix in range(Lx):
            for iy in range(Ly):
                index = ix * Ly + iy + Lx * Ly * auxT
                Ux0 = Ux[index]
                Uy0 = Uy[index]
                rho0 = self.rho(ix, iy, True)
                for i in range(Q):
                    int n0 = self.n(ix, iy, i)
                    self.fnew[n0] = self.feq(rho0, Ux0, Uy0, i)

    cpdef Advection(self):
        cdef int ix, iy, i, ixnext, iynext, n0, n0next
        for ix in range(Lx):
            for iy in range(Ly):
                for i in range(Q):
                    ixnext = ix + self.Vx[i]
                    iynext = iy + self.Vy[i]
                    if 0 <= ixnext < Lx and 0 <= iynext < Ly:
                        n0 = self.n(ix, iy, i)
                        n0next = self.n(ixnext, iynext, i)
                        self.f[n0next] = self.fnew[n0]

    cpdef PrintData(self, str NameFile, double t):
        cdef FILE* MyFile = fopen(NameFile.encode('utf-8'), "w")
        cdef double rho0, Ux0, Uy0
        cdef int step = 1
        for ix in range(0, Lx, step):
            for iy in range(0, Ly, step):
                rho0 = self.rho(ix, iy, False)
                Ux0 = self.Jx(ix, iy, False) / rho0
                Uy0 = self.Jy(ix, iy, False) / rho0
                fprintf(MyFile, "%d %d %f %f %f\n", ix, iy, rho0, Ux0, Uy0)
            fprintf(MyFile, "\n")
        fclose(MyFile)


import sys
from libc.stdio cimport printf, system
from libc.stdlib cimport malloc, free
from libcpp.string cimport stringstream

# Importa la clase LatticeBoltzman desde el archivo generado por Cython
from lattice_boltzmann cimport LatticeBoltzman, LoadData

def main():
    cdef int t, tframe = 200, tmax, delta_t = 1
    cdef double rho0 = 10.0, Ux0 = 0.0, Uy0 = 0.0
    cdef double mu_x, mu_y, sigma_x, sigma_y
    cdef double D
    cdef double tau, Utau, UmUtau
    cdef int ret

    # Parámetros generales de la simulación
    tmax = t_hour * 40

    # Parámetros para la distribución gaussiana que inicializa la densidad
    mu_x = Lx / 2.0
    mu_y = Ly / 2.0
    sigma_x = Lx / 32.0
    sigma_y = Ly / 32.0

    # Crear una instancia de la clase LatticeBoltzman
    cdef LatticeBoltzman Air = LatticeBoltzman()

    # Cargar los datos de velocidad
    LoadData("velocity.txt")

    # Leer parámetros desde la línea de comandos
    if len(sys.argv) > 1:
        D = float(sys.argv[1])  # Coeficiente de difusión
    else:
        printf("Uso: python main.py <coeficiente_de_difusion>\n")
        return

    # Calcular el tiempo de relajación tau basado en el coeficiente de difusión y el paso de tiempo
    tau = (D / delta_t * Cs2) + 0.5

    # Calcular otros valores útiles basados en tau
    Utau = 1.0 / tau
    UmUtau = 1 - Utau  # 1 - 1/tau

    # Inicializar la simulación con las condiciones iniciales (densidad y velocidades)
    Air.Start(rho0, Ux0, Uy0, mu_x, mu_y, sigma_x, sigma_y)

    # Bucle principal de la simulación
    for t in range(tmax + 1):
        Air.Collision()
        Air.ImposeFields(t)
        Air.Advection()

        # Guardar resultados cada tframe pasos
        if t % tframe == 0:
            # Mostrar el porcentaje de avance de la simulación en la consola
            printf("Porcentaje de avance: %d%%\n", (t * 100) / tmax)

# Solo ejecuta main si se llama desde la línea de comandos
if _name_ == "_main_":
    main()
    # from setuptools import setup
# from Cython.Build import cythonize
# import numpy

# setup(
#     ext_modules=cythonize("lattice_boltzmann.pyx", compiler_directives={'language_level': "3"}),
#     include_dirs=[numpy.get_include()]
# )
