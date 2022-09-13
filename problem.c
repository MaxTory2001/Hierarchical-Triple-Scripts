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
void generate_sims(struct reb_simulation* r);
void global_orbit_output(struct reb_simulation* r, struct reb_particle p);
void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p);
void make_orbit();
bool stopped;
double tmax = 2.*M_PI*1.0e1;     // 100 years translated into code units
int n = 0; float inc; float a_over_rh; float mass_ratio = 0.0;
float ecc = 0.0, a_out;
// For plotting info about a specific orbit
int target_inc = 50;
float target_a = 0.4;

// hill radius
float r_hill;

int main(int argc, char* argv[]){
    // Initial conditions
    if (argc > 1){
        // Mass ratio between planet and outer star
        mass_ratio = pow(10, strtof(argv[1], NULL));
	a_out = pow(1 + mass_ratio, 1.0/3.0);
	printf("a_out = %.3f\n", a_out);
        r_hill = pow((mass_ratio)/3., 1.0/3.0) * a_out;
	printf("r_hill = %.3f\n", r_hill);
        if (argc == 4) {
            // Only check one given inclination - for computation time
            inc = strtof(argv[2], NULL);
            a_over_rh = strtof(argv[3], NULL);
                
            make_orbit(); 
        } 
    }
    else {
        printf("Must provide mass ratio as command line argument!\n");
        return 1;
    }    
}

void make_orbit(){
    /*
     *
    */
        
        struct reb_simulation* r = reb_create_simulation();
        // Setup constants
        r->dt             = M_PI*1.0e-3;     // initial timestep
        r->integrator        = REB_INTEGRATOR_IAS15;
        r->heartbeat        = heartbeat;

        // Adding star
        struct reb_particle planet = {0}; 
        planet.m  = mass_ratio;
        planet.hash = reb_hash("m1");
	printf("m1 = %d", planet.hash);
        reb_add(r, planet); 

        generate_sims(r);

        n++;

        reb_free_simulation(r);
}

void generate_sims(struct reb_simulation* r){
    // 0 mass particle - should have no energy transfer
    for (float M = 0; M < 2 * M_PI; M += M_PI/6){
        struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0], 0, a_over_rh * r_hill, ecc, inc * M_PI / 180, 0.0, 0.0, reb_tools_M_to_f(M));
        m2.hash = reb_hash("m2");
        printf("\nm2 = %d", m2.hash);
        reb_add(r, m2);

        // The star
        struct reb_particle star = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0, a_out, 0.0, 0.0, 0.0, 0.0, 0.0);
        reb_add_fmt(r, "m a e inc", 1.0, a_out, 0.0, inc); 

        reb_move_to_com(r);
        
        printf("\ninc = %.4f, a = %.4f, n = %d\n", inc, a_over_rh, n);
        stopped = false;
        reb_integrate(r, tmax);
        // when orbit was either finished or disrupted, output time
    }
}

void heartbeat(struct reb_simulation* r){
    float d_in = reb_particle_distance(&r->particles[0], &r->particles[1]);
    float d_out = reb_particle_distance(&r->particles[0], &r->particles[2]);
    // diagnostics about 1 in every 199 orbits

    if(reb_output_check(r, 0.01)){        // outputs to the screen
    	char orbits[64];

    	sprintf(orbits, "orbits/orbits_all_M_%.2e_i_%.1f.txt", mass_ratio, inc);

    	reb_output_orbits(r, orbits);

    	relative_orbit_output(r, r->particles[1], r->particles[2]);
    	global_orbit_output(r, r->particles[0]);
    	global_orbit_output(r, r->particles[1]);
    	global_orbit_output(r, r->particles[2]);
    }
    
    // exiting when each orbit is disrupted (more than 3* the start energy or energy has somehow gone negative)
    if(d_out/2 < d_in && !stopped){          
        stopped = true;
	print_time(r);
    }
}

void print_time(struct reb_simulation* r){
    char fname[64];

    sprintf(fname, "orbits/time_all_M_%.2e_i_%.1f_a_%.2f.txt", mass_ratio, inc, a_over_rh);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.7f\t%d\t%d\n", r->t, n, (int) inc);
    fclose(fp);
}

void global_orbit_output(struct reb_simulation* r, struct reb_particle p) {
    // prints detailed info about an orbit - orbital elements and more
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, p, r->particles[2]);
    char fname[64];
    sprintf(fname, "orbits/detailed_global_%.2e_%.2f_%.2f_%3d.txt", mass_ratio, inc, a_over_rh, p.hash);
    FILE* fp = fopen(fname, "a");
    fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, p.x, p.y, p.z, p.vx, p.vy, p.vz, o.a, o.e, o.inc, o.Omega);
    fclose(fp);
}

void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p) {
    /* 
    * p: primary particle - centre of the orbit
    * s: secondary particle
    * prints detailed info about an orbit - orbital elements and more
    **/
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, s, p);
    char fname[64];
    sprintf(fname, "orbits/detailed_relative_%.2e_%.2f_%.2f.txt",mass_ratio, inc, a_over_rh); 
    FILE* fp = fopen(fname, "a");
    fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, (s.x-p.x) / r_hill, (s.y-p.y) / r_hill, (s.z-p.z) / r_hill, s.vx-p.vx, s.vy-p.vy, s.vz-p.vz, o.a, o.e, o.inc, o.omega);
}
