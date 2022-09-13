"""
Max Tory
Evolution of one orbit which is stable by the Tory-Grishin-Mandel criterion
and unstable by the Vynatheya-Hamers criterion. Orbital elements are tracked
to follow the evolution of this orbit and when it becomes unstable by each
criterion
"""

import rebound
import sys
import math

class System:
    def __init__(self):
        self.sim = None
        self.e_in, self.e_out = 0.0, 0.0
        self.previous_output = -1
        self.theta = 0

        if len(sys.argv) == 5:
            argv = sys.argv
            self.q_out = 10 ** float(argv[1])
            self.q_in = 10 ** float(argv[4])
            self.r_h = (self.q_out/3) ** (1./3.)
            self.a_out = (1 + self.q_out) ** (1./3.)
            self.inc = float(argv[2])
            self.a_over_rh = float(argv[3])
        else:
            print("Wrong number of args")

    def make_orbit(self, theta):
        self.sim.add(m=self.q_out/(1 + self.q_in))
        self.sim.add(m=self.sim.particles[0].m * self.q_in,
                a=self.a_over_rh * self.r_h,
                e=self.e_in,
                inc=self.inc * math.pi/180,
                theta=theta)
        self.sim.add(m=1,
                a=self.a_out,
                e=self.e_out)
        self.sim.move_to_com()
        self.sim.integrate(200 * math.pi)

    def setup(self):
        for theta in [math.pi/10 * i for i in range(20)]:
            self.sim = rebound.Simulation()
            self.sim.G = 1.0;
            self.sim.dt = math.pi * 1e-10
            self.sim.heartbeat = self.heartbeat

            self.theta=theta
            try:
                self.make_orbit(self.theta)
            except e:
                pass
            self.sim = None
            self.previous_output = -1


    def heartbeat(self, sim_pointer):
        if self.sim.t >= self.previous_output + 1:
            orbits = self.sim.calculate_orbits()
            self.previous_output = self.sim.t
            with open(f"track_orbit_f_{self.theta:,.2f}.txt", "a") as f:
                f.writelines([f"{self.sim.t:.2f}\t{o.a:.4f}\t{o.e:.4f}\t{o.inc:.2f}\t{o.f:.4f}\n" for o in orbits])


def main():
    system = System()
    system.setup()


if __name__ == "__main__":
    main()

