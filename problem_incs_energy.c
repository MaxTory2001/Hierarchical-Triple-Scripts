/**
 * Hill Stability Breakdown
 * 
 * Here I use the IAS15 integrator to check stability of orbits
 * in heirarchical triples for different inclinations and semi-major axes,
 * then extend this to less hierarchical regimes to determine the breakdown
 * point of the hill stability criterion
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "rebound.h"
#include <time.h>

void heartbeat(struct reb_simulation* r);
void print_time(struct reb_simulation* r);
void generate_sims(struct reb_simulation* r, double r1, double l1);
void global_orbit_output(struct reb_simulation* r, struct reb_particle p);
void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p);
void loop_sma();

float a_out; double a_over_aout;
double tmax = 2.*M_PI*1.0e2;     // 100 years translated into code units
int n = 0; float inc; float dE; float mass_ratio = 0.0; float stop_inc = 0.0;
float ecc = 0.0;
double f = 0.0;
// hill radius
float r_hill;

int main(int argc, char* argv[]){
    // Initial conditions
    if (argc > 1){
        // Mass ratio between planet and outer star
        mass_ratio = pow(10, strtof(argv[1], NULL));

        r_hill = pow((mass_ratio)/3., 1.0/3.0);
	a_out = pow(1 + mass_ratio, 1.0/3.0);
        if (argc == 3) {
            // Only check one given inclination - for computation time
            inc = strtof(argv[2], NULL);
            stop_inc = inc + 5.0;
            // Seeding random number generator
	    srandom((int) inc);
            while (inc < stop_inc) {
                loop_sma();
                inc ++;
            }
                
        } else {
            for (f = 0.0; f <= 180.0; f++){
                loop_sma();
            }
        }
    } else {
        printf("Must provide mass ratio as command line argument!\n");
        return 1;
    }
}

// ~~~~~~~~~~~~~~~~~ Functions for root finding magic ~~~~~~~~~~~~~~~~~~

double mag(struct reb_vec3d vec) {
    return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

double effective_potential(struct reb_vec3d p1, struct reb_vec3d p2, struct reb_vec3d R_vec) {
    // p1 is position of planet, p2 is position of star
    double R = mag(R_vec);
    struct reb_vec3d r1_vec = {.x = R_vec.x - p1.x,
                                .y = R_vec.y - p1.y,
                                .z = R_vec.z - p1.z};
    struct reb_vec3d r2_vec = {.x = R_vec.x - p2.x,
                                .y = R_vec.y - p2.y,
                                .z = R_vec.z - p2.z};
    double r1 = mag(r1_vec), r2 = mag(r2_vec);
    return - mass_ratio/r1 - 1/r2 - 0.5*R*R;
}

double effective_potential_grad_grad(double m1, double m2, double x){
    double a1 = (m2/(m1+m2))*pow(m1 + m2, 1./3.);
    double a2 = -(m1/(m1+m2))*pow(m1 + m2, 1./3.);

    return -2*m1/fabsl(pow((x - a1), 3)) - 2*m2/fabsl(pow((x - a2), 3)) - 1;
}

double effective_potential_grad(double m1, double m2, double x){
    double a1 = (m2/(m1+m2))*pow(m1 + m2, 1./3.);
    double a2 = -(m1/(m1+m2))*pow(m1 + m2, 1./3.);
    return m1/(fabsl((x-a1)) * (x-a1)) + m2/(fabsl((x-a2)) * (x-a2)) - x;
}

double total_minus_L1(struct reb_vec3d p1, struct reb_vec3d p2, double r1, double U_L1) {
    struct reb_vec3d R_vec;
    R_vec.x = p1.x + r1*cos(f);
    R_vec.y = p1.y + r1*sin(f)*cos(inc);
    R_vec.z = p1.z + r1*sin(f)*sin(inc);

    return effective_potential(p1, p2, R_vec) + 0.5*mass_ratio/r1 - U_L1;
}

double find_L1(double m1, double m2, double x_old, double xmin, double xmax, double dU_dx_old, double err) {
    // Finding the L1 lagrange point using Newton's method to find the root of the 
    // gradient of the effective potential function along the x axis
    if (m1 == m2) return 0;
    else if (m1 < m2) printf("m1 cannot be < m2");
    double x = x_old - dU_dx_old/effective_potential_grad_grad(m1, m2, x_old);
    // should never happen but just in case we are out of bounds, guess back in the middle
    if (x <= xmin || x >= xmax) return find_L1(m1, m2, 0.5*(xmin + xmax), xmin, xmax, 0, err);
    double dU_dx = effective_potential_grad(m1, m2, x);

    if (fabsl(dU_dx) < err) return x;
    else return find_L1(m1, m2, x, xmin, xmax, dU_dx, err);
}

double find_root_r1(struct reb_vec3d p1, struct reb_vec3d p2, double lower, double upper, double threshold, double U_eff, int count) {
    // finding the r value that will give a desired total energy using the bisection method
    if (upper < lower || count > 50) {
	printf("Didn't converge! No root in range!");
	return -50000;
    }
    double middle = 0.5 * (lower + upper);
    double e1 = total_minus_L1(p1, p2, lower, U_eff);
    double e2 = total_minus_L1(p1, p2, upper, U_eff);
    if (e1 * e2 > 0) {
	return find_root_r1(p1, p2, lower, upper-0.1, threshold, U_eff, count + 1);
    }
    double e3 = total_minus_L1(p1, p2, middle, U_eff);
    if (fabsl(e3) <= threshold) {
	return middle;
    }
    else if (e3/e1 > 0) {
	return find_root_r1(p1, p2, middle, upper, threshold, U_eff, count + 1);
    }
    
    else{
	return find_root_r1(p1, p2, lower, middle, threshold, U_eff, count + 1);
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void loop_sma(){
    /*
     *
    */
    // getting lagrange point
    double l1 = find_L1(1.0, mass_ratio, 0, -1/(1 + mass_ratio), 0, 0, 1e-6);
    // getting energy at lagrange point
    struct reb_vec3d pos = { .x = l1, .y = 0, .z = 0};
    struct reb_vec3d p1 = { .x = -1/(1 + mass_ratio)*a_out, .y = 0, .z = 0};
    struct reb_vec3d p2 = { .x = mass_ratio/(1 + mass_ratio)*a_out, .y = 0, .z = 0};
    double U_eff = effective_potential(p1, p2, pos);
    for (dE = -0.1 * U_eff; dE <= 0.6 * U_eff; dE += 0.005){
        
        struct reb_simulation* r = reb_create_simulation();
        // Setup constants
        r->dt             = M_PI*1.0e-3;     // initial timestep
        r->integrator        = REB_INTEGRATOR_IAS15;
        r->heartbeat        = heartbeat;

        // Adding star
        struct reb_particle planet = {0}; 
        planet.m  = mass_ratio;
        planet.hash = reb_hash("m1");
        reb_add(r, planet); 

        // random mean anomoly
	//f = 2 * M_PI * ((double)random())/RAND_MAX;
        // getting radius of the particle from the planet to have the desired energy
        double r1 = find_root_r1(p1, p2, 1e-7, r_hill, 1e-4, U_eff + dE, 0);
	if (r1 < 0) {
	    print_time(r);
	    n++;
	    reb_free_simulation(r);
	    continue;
	}

        generate_sims(r, r1, l1);

        // 1 in ~100 orbits have their diagnostics saved. 97 is prime so we get a good combination of inclinations and a's
        if (n % 199 == 0){
            char fname[64];

            if (stop_inc == 0.0){
                sprintf(fname, "energy_incs/diagnostics_0_M%.2e.txt", mass_ratio);
            } else {
                sprintf(fname, "energy_incs/diagnostics_0_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);
            }

            FILE* fp = fopen(fname, "a");
            fprintf(fp, "\n");
            fclose(fp);
        }

        n++;

        reb_free_simulation(r);
    }
}

void generate_sims(struct reb_simulation* r, double r1, double l1){
    a_over_aout = r1/a_out;
    printf("r1 = %.4f\n", r1);
    printf("a_out = %.4f", a_out);
    // the satelite (m2)
    struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0], 0.0, r1, ecc, inc * M_PI / 180, 0.0, 0.0, f * M_PI / 180);
    m2.hash = reb_hash("m2");
    reb_add(r, m2);

    // The star
    struct reb_particle star = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0, a_out, 0.0, 0.0, 0.0, 0.0, 0.0);
    star.hash = reb_hash("star"); 

    reb_add(r, star); 

    reb_move_to_com(r);

    printf("\nf = %.4f, dE = %.4f, a = %.4f, n = %d, l = %.4f\n", f, dE, a_over_aout / r_hill, n, l1);

    char orbits[64];

    if (stop_inc == 0.0){
        sprintf(orbits, "energy_fs/orbits_0_M%.2e.txt", mass_ratio);
    } else {
        sprintf(orbits, "energy_fs/orbits_0_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);
    }

    reb_output_orbits(r, orbits);

    reb_integrate(r, tmax);
    // when orbit was either finished or disrupted, output time
    print_time(r);
}

void heartbeat(struct reb_simulation* r){
    float d_in = reb_particle_distance(&r->particles[0], &r->particles[1]);
    float d_out = reb_particle_distance(&r->particles[1], &r->particles[2]);
    // diagnostics about 1 in every 199 orbits
    if(reb_output_check(r, 10*M_PI)){        // outputs to the screen
        if (n % 199 == 0){
            char fname[64];

            if (stop_inc == 0) {
                sprintf(fname, "energy_fs/diagnostics_0_M%.2e.txt", mass_ratio);
            } else {
                sprintf(fname, "energy_fs/diagnostics_0_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);

            }
            
            FILE* fp = fopen(fname, "a");
            struct reb_vec3d ang_mom = reb_tools_angular_momentum(r);
            fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, reb_tools_energy(r), ang_mom.x, ang_mom.y, ang_mom.z);
            fclose(fp);
        }
    }
    
    // exiting when each orbit is disrupted (more than 3* the start energy or energy has somehow gone negative)
    if(d_out/2 < d_in){          
        r->status = 1;
    }
}

void print_time(struct reb_simulation* r){
    char fname[64];

    sprintf(fname, "energy_fs/time_0_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.7f\t%d\t%d\t%d\n", r->t, n, (int) f, (int) stop_inc);
    fclose(fp);
}
/*
void global_orbit_output(struct reb_simulation* r, struct reb_particle p) {
    // prints detailed info about an orbit - orbital elements and more
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, p, r->particles[2]);
    FILE* fp = fopen("detailed_global.txt", "a");
    fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, p.x, p.y, p.z, p.vx, p.vy, p.vz, o.a, o.e, o.inc, o.Omega);
    fclose(fp);
}

void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p) {
    
    * p: primary particle - centre of the orbit
    * s: secondary particle
    
    prints detailed info about an orbit - orbital elements and more
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, s, p);
    FILE* fp = fopen("detailed_relative.txt", "a");
    fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, (s.x-p.x) / r_hill, (s.y-p.y) / r_hill, (s.z-p.z) / r_hill, s.vx-p.vx, s.vy-p.vy, s.vz-p.vz, o.a, o.e, o.inc, o.omega);
    fclose(fp);
}*/
