# -*- coding: utf-8 -*-
"""Ejercicio2-d

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1bcthC3uLPGMkPZk4Hw2_a-Uqv_gS83kq
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import integrate

#Factor de amplificación
#La gráficas que vienen con el código usaron FA=10000
FA=100

#La solución general de la ecuación diferencial planteada en el ejercicio 2 es la función de Bessel para n=0 de λ*r
#R_λ(r)=J_0(λ*r)
def Bessel(n,l,r):
  t=np.linspace(0, np.pi, FA)
  return (1/np.pi)*integrate.simpson(np.cos(n*t-l*r*np.sin(t)), t)

#l: Arreglo con la información de λ
l=np.linspace(0, 15.0, 15*FA)
#f_l: Arreglo con la información de f(λ)
f_l=np.zeros(len(l))
#f(λ)=R(r=1;λ)=J_0(λ)
for i in range(len(l)):
  f_l[i]=Bessel(0, l[i], 1)

#Creamos y guardamos la gráfica de (λ, f(λ))
fig = plt.figure(figsize=(17, 10.5))
plt.title('Solucion teorica de f(λ)', fontsize=20)
plt.xlabel('$λ$', fontsize=15)
plt.ylabel('$f(λ)=J_0(λ)$', fontsize=15)
plt.plot(l, f_l)
plt.savefig('Grafica_2d1.png')
plt.close()

#Esta función nos devuelve el i-ésimo término del arreglo l asociado a un valor del dominio λ ([0,15])
def i(x):
  return int(x*FA)

#Con dos valores del dominio A, B, esta función halla el valor M para el cual f(M) se hace cero con cierto error
def CeroPorBiseccion(data, A, B, Err):
  a=i(A)
  b=i(B)
  while B-A>Err:
    M=(A+B)/2
    m=i(M)
    if (data[a]*data[m])>0:
      A=M
    else:
      B=M
  return M

#Esta función grafica la función de Bessel para un determinado λ_i
def Graficar_integral_Bessel(Lambda_i):
  #Dominio: [0, 1.0]
  r=np.linspace(0, 1.0, 10*FA)

  #Función de Bessel evaluada en cada punto del dominio (rango)
  Bessel_plot=np.zeros(len(r))
  for i in range(len(r)):
    Bessel_plot[i]=Bessel(0, Lambda_i, r[i])

  #Graficamos J_0(λ*r) con una leyenda que contiene el valor de λ usado acotado a tres decimales
  plt.plot(r, Bessel_plot, label=(u"λ={}".format(round(Lambda_i,3))))
  plt.legend(fontsize="16")

#Esta función grafica los modos normales que no dependen de la coordenada angular 𝜃, para ello
#halla los ceros de f(λ) cuya solución teórica fue hallada anteriormente y los usa como los λ_i que
#tienen asociados cada uno la función de Bessel J_0(λ_i*r), y superpone las gráficas
def Graficador_modos_normales(data_x, data_y, N, Err):
  #Arreglo que contiene los λ_i
  Lambda=np.zeros(N-1)

  #Arreglo que contiene los puntos que se obtienen de particionar el dominio en N-1 partes
  #(cada parte para cada cero a hallar en f(λ))
  Particion=np.linspace(data_x[0], data_x[-1], N)

  #Título de la gráfica y sus ejes, y su tamaño
  fig = plt.figure(figsize=(17, 10.5))
  plt.title('MODOS NORMALES (BESSEL)', fontsize=20)
  plt.xlabel('$r$', fontsize=15)
  plt.ylabel('$J_0(λ \cdot r)$', fontsize=15)

  #En este loop hallamos cada λ_i yendo por cada una de las particiones del dominio y graficamos su solución asociada
  for i in range(N-1):
    Lambda[i]=CeroPorBiseccion(data_y, Particion[i], Particion[i+1], Err)
    Graficar_integral_Bessel(Lambda[i])

  #Guardamos y cerramos la gráfica
  plt.savefig('Grafica_2d2.png')
  plt.close()

#Usamos la función anterior para graficar todos los modos normales
Graficador_modos_normales(l, f_l, 6, 0.001)
