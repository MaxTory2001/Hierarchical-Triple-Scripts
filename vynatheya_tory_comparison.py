"""
Max Tory
Stability of a map of orbits against a_in and inc according to different criteria
"""

import rebound
import sys
import math
import numpy as np

class System:
    def __init__(self):
        self.sim = None
        self.e_in, self.e_out = 0.0, 0.0
        self.previous_output = -1
        self.theta = 0

        if len(sys.argv) == 4:
            argv = sys.argv
            self.q_out = 10 ** float(argv[1])
            self.q_in = 10 ** float(argv[3])
            self.r_h = (self.q_out/3) ** (1./3.)
            self.a_out = (1 + self.q_out) ** (1./3.)
            self.inc = float(argv[2])
            self.stop_inc = self.inc + 5
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

    def setup_systems(self):
        """Setting up all the orbits at different inclinations and a_ins"""
        while self.inc < self.stop_inc:
            for self.a_over_rh in np.arange(0.2, 1.5, 0.01):
                stable = 0
                hamers_stable = 0
                num_angles = 20
                for theta in [math.pi * 2/num_angles * i for i in range(num_angles)]:
                    self.sim = rebound.Simulation()
                    self.sim.G = 1.0;
                    self.sim.dt = math.pi * 1e-10
                    self.sim.heartbeat = self.heartbeat

                    self.theta=theta
                    try:
                        self.make_orbit(self.theta)
                    except e:
                        pass
                    if self.sim.status != 1:
                        stable += 1
                        if self.vynatheya_hamers_criterion(self.sim):
                            hamers_stable += 1

                    self.sim = None
                    self.previous_output = -1

                self.print_stable(stable/num_angles, hamers_stable/num_angles)
            self.inc += 1

    def print_stable(self, stable, hamers_stable):
        path = f"hamers/stable_orbits_q_{self.q_out:,.2e}_i_{self.inc:,.1f}_q_in_{self.q_in:,.2e}.txt"
        hamers_path = f"hamers/hamers_stable_orbits_q_{self.q_out:,.2e}_i_{self.inc:,.1f}_q_in_{self.q_in:,.2e}.txt"
        with open(path, "a") as f:
            f.write(f"{stable:.2f}\n")

        with open(hamers_path, "a") as f:
            f.write(f"{hamers_stable:.2f}\n")

    def heartbeat(self, sim_pointer):
        p1, p2 = self.sim.particles[:2]
        d_in = p1 ** p2
        if (d_in > self.a_out / 2):
            self.sim.status = 1


    def vynatheya_hamers_criterion(self, sim):
        orbits = sim.calculate_orbits()
        final_a_in, final_a_out = orbits[0].a, orbits[1].a

        a_in = self.a_over_rh * self.r_h * self.a_out
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
    
