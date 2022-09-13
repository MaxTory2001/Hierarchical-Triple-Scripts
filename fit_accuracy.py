"""
Max Tory
Stability of a random selection of orbits against a_in and inc according to
different criteria. For evaluating fits to the stability boundary
"""

import rebound
import sys
import math
import random
import numpy as np

class System:
    def __init__(self):
        self.sim = None
        self.e_in, self.e_out = 0.0, 0.0
        self.previous_output = 0
        self.theta = 0

        if len(sys.argv) == 1:
            self.seed = 1
            random.seed(1)
        else:
            self.seed = int(sys.argv[1])
            random.seed(int(sys.argv[1]))

    def make_orbit(self):
        self.sim.add(m=self.q/(1 + self.q_in))
        self.sim.add(m=self.sim.particles[0].m * self.q_in,
                a=self.a_over_rp * (1 - self.e_out) * self.a_out,
                e=self.e_in,
                inc=self.inc,
                M=self.M1,
                Omega=self.Omega,
                omega=self.omega)
        self.sim.add(m=1,
                a=self.a_out,
                e=self.e_out,
                M=self.M2)
        self.sim.move_to_com()
        print("integrating")
        try:
            while self.sim.t < 200 * math.pi and self.sim._status != 1:
                self.sim.integrate(0.1 * math.pi)
                self.heartbeat()
        except:
            pass

    def predict_stable_a_over_rp(self):
        f = 10 ** (-0.52 + 0.04 * self.q) * self.q ** (0.32 + 0.1 * self.q)
        a0, a1, a2, a3, a4 = -0.0798, 0.5045, -0.8617, 0.2260, 0.7300
        g_poly = a0 * pow(self.inc, 4) + a1 * pow(self.inc, 3) + a2 * pow(self.inc, 2) + a3 * self.inc + a4
        g = g_poly if self.inc >= math.pi/3 else -0.17 * np.cos(self.inc) + 0.63 
        g /= 0.54
        h = 10 ** ((self.inc * 180/math.pi) * (self.q ** 1.3)/1500)

        return f * g / h

    def setup_systems(self):
        """Setting up all the orbits at different inclinations and a_ins"""
        for i in range(4000):
            # vynatheya test

            self.sim = rebound.Simulation()
            self.sim.G = 1.0;
            self.sim.dt = 1e-10
            self.sim.heartbeat = self.heartbeat

            self.inc = math.acos(random.random() * 2 - 1)
            self.e_in, self.e_out = random.random(), random.random()
            self.M1, self.M2 = 2 * math.pi * random.random(), 2 * math.pi * random.random()
            self.Omega, self.omega = 2 * math.pi * random.random(), 2 * math.pi * random.random()

            self.q = 10 ** (-6 * random.random())
            self.q_in = 10 ** (-2 * random.random())

            self.r_h = (self.q / 3) ** (1./3.)
            self.a_over_rp = (random.random() + 0.5) * self.predict_stable_a_over_rp()
            self.a_out = (1 + self.q) ** (1./3.)

            print(f'q = {self.q}, inc={self.inc}, a={self.a_over_rp/self.predict_stable_a_over_rp()}, e_in={self.e_in}, e_out={self.e_out}')

            self.make_orbit()
            with open(f"fit_accuracy/stability_predictions_seed_{self.seed}.txt", "a") as f:
                f.write(f"{math.log10(self.q):.3f}\t{math.log10(self.q_in):.3f}\t{(self.inc * 180 / math.pi):.5f}\t{self.e_in:.5f}\t{self.e_out:.5f}\t{self.a_over_rp:.5e}\t{(self.sim._status + 1) % 2}\n")
            self.sim = None


    def heartbeat(self, sim_pointer=None):

        p1, p2 = self.sim.particles[:2]
        d_in = p1 ** p2

        if (d_in > self.a_out * (1-self.e_out) / 2):
            self.sim._status = 1
            print("unstable")

    def vynatheya_hamers_criterion(self, sim):
        orbits = sim.calculate_orbits()
        final_a_in, final_a_out = orbits[0].a, orbits[1].a

        a_in = self.a_over_rp * (1 - self.e_out) * self.a_out
        a_out = self.a_out

        stable = True
        if abs((a_in - final_a_in)/a_in) > 0.1: stable = False
        elif abs((a_out - final_a_out)/a_out) > 0.1: stable = False
        return stable


def main():
    system = System()
    system.setup_systems()


if __name__ == "__main__":
    main()
