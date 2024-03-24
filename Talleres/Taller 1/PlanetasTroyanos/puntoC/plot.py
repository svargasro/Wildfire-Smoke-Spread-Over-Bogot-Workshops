import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


x_s, y_s, x_j, y_j, x_t, y_t, t = np.genfromtxt('data.txt', unpack=True, usecols=(0, 1, 2, 3, 4, 5, 6))



plt.style.use('seaborn-v0_8')

with PdfPages('plot.pdf') as pdf:

    fig, axes = plt.subplots(4,1,figsize=(7, 7))

    axes[0].plot(t, x_s, '.', color='yellow', label=r"$x^{'}(t)_{S}$")
    axes[1].plot(t, y_s, '.', color='yellow', label=r"$y^{'}(t)_{S}$")
    axes[2].plot(t, x_j, '.', color='black', label=r"$x^{'}(t)_{J}$")
    axes[3].plot(t, y_j, '.', color='black', label=r"$y^{'}(t)_{J}$")


    # Se ajustan demás detalles del gráfico.
    for i in [0,2]:
        axes[i].set_xlabel('t', fontsize=12)
        axes[i].set_ylabel('x\'(t)', fontsize=12)
        axes[i].legend(loc='upper left')
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
    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


    fig, axes = plt.subplots(3,1,figsize=(7, 7))

    axes[0].plot(t, x_t, '.', color='cyan', label=r"$x^{'}(t)_{T}$")
    axes[1].plot(t, y_t, '.', color='cyan', label=r"$y^{'}(t)_{T}$")


    axes[0].set_xlabel('t', fontsize=12)
    axes[0].set_ylabel('x\'(t)', fontsize=12)
    axes[0].legend(loc='upper right')
    axes[0].grid(True, linestyle='--')
    axes[0].set_title("x' vs. t'", fontsize=14)
    axes[0].ticklabel_format(style='sci', axis='x', scilimits=(0,0))


    axes[1].set_xlabel('t', fontsize=12)
    axes[1].set_ylabel('y\'(t)', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, linestyle='--')
    axes[1].set_title("y' vs. t", fontsize=14)
    axes[1].ticklabel_format(style='sci', axis='x', scilimits=(0,0))



    axes[2].plot(x_t, y_t, ".", color="cyan", label=r"$(x',y')_{T}$")
    axes[2].set_xlabel("x'", fontsize=12)
    axes[2].set_ylabel("y'", fontsize=12)
    axes[2].legend(loc='upper right')
    axes[2].grid(True, linestyle='--')
    axes[2].set_title("y' vs. x'", fontsize=14)
    axes[2].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
