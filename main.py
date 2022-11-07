import matplotlib.pyplot as plt
import numpy as np

from matplotlib.animation import FuncAnimation


class Phys:
    def __init__(self):
        self.c = 3e8
        self.m_0 = 9.10938e-31
        self.q = 1.60217663e-19
        self.eps_0 = 8.854187e-12
        self.E = 1
        self.alpha = 45
        self.v = 1
        self.p_0 = self.m_0 / np.sqrt(1 - (self.v / self.c) ** 2) * self.v

        #  modeling consts
        self.time_limit = 1
        self.time = None
        self.x = None
        self.y = None
        self.dx_dt = None
        self.dy_dt = None

    def prepare_data(self):
        self.time = np.linspace(0, self.time_limit, 100)

        self.x = lambda t: (self.p_0 * (self.c ** 2) * np.sin(self.alpha)) / (self.q * self.E) * np.arcsinh(
            (self.c * self.q * self.E * t) / self.eps_0)

        self.y = lambda t: 1 / (self.q * self.E) * np.sqrt(
            self.eps_0 ** 2 + (self.c * self.q * self.E * t) ** 2) - self.eps_0 / (
                                   self.q * self.E) + self.p_0 * (self.c ** 2) * np.cos(self.alpha) / (
                                       self.q * self.E) * np.arcsinh(self.c * self.q * self.E * t / self.eps_0)

        self.dx_dt = lambda t: (self.p_0 * (self.c**2) * np.sin(self.alpha)) / np.sqrt(self.eps_0 ** 2 + (self.c * self.q * self.E * t) ** 2)

        self.dy_dt = lambda t: ((self.c ** 2) * (self.q * self.E * t + self.p_0 * np.cos(self.alpha))) / np.sqrt(self.eps_0 ** 2 + (self.c * self.q * self.E * t) ** 2)

    def plot(self):
        plt.grid(True)
        plt.xlabel('Время, c')
        plt.ylabel('Координата X, м')
        plt.title(r'График зависимости $X(t)$')
        plt.plot(self.time, self.x(self.time), lw=3)
        plt.savefig(f'x_coord_{self.alpha}.png', dpi=700)

        plt.close()

        plt.grid(True)
        plt.xlabel('Время, c')
        plt.ylabel('Координата Y, м')
        plt.title(r'График зависимости $Y(t)$')
        plt.plot(self.time, self.y(self.time), lw=3)
        plt.savefig(f'y_coord_{self.alpha}.png', dpi=700)

        plt.close()

        plt.grid(True)
        plt.xlabel('Время, с')
        plt.ylabel(r'Модуль скорости, $\frac{м}{с}$')
        plt.title(r'График зависимости модуля $V(t)$')
        plt.plot(self.time, np.sqrt((self.dx_dt(self.time) ** 2) + (self.dy_dt(self.time) ** 2)), lw=3)
        plt.savefig(f'v_t_{self.alpha}')


mover = Phys()
mover.prepare_data()
mover.plot()
