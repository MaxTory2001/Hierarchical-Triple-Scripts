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
void loop_sma();

float a_out; double M; double f;
double tmax = 2.*M_PI*1.0e2;     // 100 years translated into code units
int n = 0; float inc; float a_over_rh; float mass_ratio = 0.0; float stop_inc = 0.0;
float ecc = 0.0;
double UNSTABLE_MEG = 5;
// hill radius
float r_hill;

int main(int argc, char* argv[]){
    // Initial conditions
    if (argc > 1){
        // Mass ratio between planet and outer star
        mass_ratio = pow(10, strtof(argv[1], NULL));

	printf("mass ratio is %.2e", mass_ratio);
        r_hill = pow((mass_ratio)/3., 1.0/3.0);
	a_out = 1.0 + pow(mass_ratio, 1.0/3.0);
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
            for (inc = 0.0; inc <= 180.0; inc++){
                loop_sma();
            }
        }
    } else {
        printf("Must provide mass ratio as command line argument!\n");
        return 1;
    }

    
}

void loop_sma(){
    /*
     *
    */
    for (a_over_rh = 0.15; a_over_rh < 1.5; a_over_rh += 0.01){
	double M = 2 * M_PI * (double) random()/ (double) RAND_MAX;

        f = reb_tools_M_to_f(0, M);
    
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
    
        generate_sims(r);
    
        // 1 in ~100 orbits have their diagnostics saved. 97 is prime so we get a good combination of inclinations and a's
        if (n % 199 == 0){
    	char fname[64];
    
    	if (stop_inc == 0.0){
    	    sprintf(fname, "incs_first_orbit/diagnostics_rand_M%.2e.txt", mass_ratio);
    	} else {
    	    sprintf(fname, "incs_first_orbit/diagnostics_rand_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);
    	}
    
    	FILE* fp = fopen(fname, "a");
    	fprintf(fp, "\n");
    	fclose(fp);
        }
    
        n++;
    
        reb_free_simulation(r);
    }
}

void generate_sims(struct reb_simulation* r){
    // random mean anomoly
    //f = 2 * M_PI * ((double) random())/ RAND_MAX;
    // the satelite (m2)
    
    struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0e-2 * r->particles[0].m, a_out * a_over_rh * r_hill, ecc, inc * M_PI / 180, 0.0, 0.0, f);
    m2.hash = reb_hash("m2");
    reb_add(r, m2);

    // The star
    struct reb_particle star = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0, a_out, 0.0, 0.0, 0.0, 0.0, 0.0);
    star.hash = reb_hash("star"); 

    reb_add(r, star); 

    reb_move_to_com(r);
    
    // add megno particles for stability analysis
    printf("\ninc = %.4f, a = %.4f, n = %d\n", inc, a_over_rh, n);

    reb_integrate(r, tmax);
    // when orbit was either finished or disrupted, output time
    print_time(r);
}

void heartbeat(struct reb_simulation* r){
    if(r->t < 3 * M_PI && reb_output_check(r, M_PI/10)){        // outputs to the screen
        /*if (n % 199 == 0){
            char fname[64];

            if (stop_inc == 0) {
                sprintf(fname, "incs_first_orbit/diagnostics_rand_M%.2e.txt", mass_ratio);
            } else 
                sprintf(fname, "incs_first_orbit/diagnostics_rand_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);

            }
            

            FILE* fp = fopen(fname, "a");
            struct reb_vec3d ang_mom = reb_tools_angular_momentum(r);
            fprintf(fp, "%.7f\t%.7e\t%.7e\t%.7e\t%.7e\n", r->t, reb_tools_energy(r), ang_mom.x, ang_mom.y, ang_mom.z);
            fclose(fp);
        }*/
	char orbits[64];

	if (stop_inc == 0.0){
	    sprintf(orbits, "incs_first_orbit/orbits_rand_M%.2e.txt", mass_ratio);
	} else {
	    sprintf(orbits, "incs_first_orbit/orbits_rand_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);
	}
	
	FILE* fp = fopen(orbits, "a");
	if (r->t == 0) fprintf(fp, "\n");
	fclose(fp);

	reb_output_orbits(r, orbits);
	
    }
    double d1 = reb_particle_distance(&r->particles[0], &r->particles[1]);

    if (d1 > a_out/1.5) {
        r->status = 1;
    }
}

void print_time(struct reb_simulation* r){
    char fname[64];

    sprintf(fname, "incs_first_orbit/time_rand_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.7f\t%d\t%d\t%d\n", r->t, n, (int) inc, (int) stop_inc);
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
