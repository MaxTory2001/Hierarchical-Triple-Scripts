"""
Script for animating a very small mass particle in orbit around a more massive
binary in a co-rotating reference frame.
"""
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from matplotlib.animation import PillowWriter

incs = (0, 30, 60 , 90, 120)
a_s = (0.2, 0.4, 0.6, 0.8, 1.0)
axes = ((1, 2), (1, 3), (2, 3))
axes_names = ("xy", "xz", "yz")

for inc in incs:
    for a_over_rh in a_s:
        for i, axis in enumerate(axes):

            fig, ax = plt.subplots(figsize=(5, 5))
            ax = plt.axes(xlim=(-1, 1), ylim=(-1, 1))
            m1, = ax.plot([], [], 'g.', markersize=20, label="m1")
            m2, = ax.plot([], [], 'r.', markersize=10, label="m2")
            m3, = ax.plot([], [], 'b.', markersize=20, label="m3")

            ax.plot(0, 0, 'X', markersize=5, color="yellow")
            plt.grid(True, lw=0.3)
            plt.legend()
            plt.title(r"$a/r_H = $%.2f, inc = %.2f, plot of %s" %
                      (a_over_rh, inc, axes_names[i]))
            hashes = [1804614826, 1973855615, -918039340]

            data_m2 = np.loadtxt(r"orbits/detailed_global_%.2f_%.2f_%d.txt" %
                                 (inc, a_over_rh, hashes[0]))
            data_m1 = np.loadtxt(r"orbits/detailed_global_%.2f_%.2f_%d.txt" %
                                 (inc, a_over_rh, hashes[1]))
            data_m3 = np.loadtxt(r"orbits/detailed_global_%.2f_%.2f_%d.txt" %
                                 (inc, a_over_rh, hashes[2]))

            def animate(i):
                m1.set_data(data_m1[i, axis[0]], data_m1[i, axis[1]])
                m2.set_data(data_m2[i, axis[0]], data_m2[i, axis[1]])
                m3.set_data(data_m3[i, axis[0]], data_m3[i, axis[1]])
                return m1, m2, m3

            anim = FuncAnimation(fig, animate, frames=6280, interval=100,
                                 repeat=False)
            anim.save(r'orbit_a=%.2f_inc=%.2f_%s.gif' %
                      (a_over_rh, inc, axes_names[i]), writer='pillow')
            plt.show()
