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
void print_megno(struct reb_simulation* r);
void generate_sims(struct reb_simulation* r);
void global_orbit_output(struct reb_simulation* r, struct reb_particle p);
void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p);
void loop_sma();
void print_stable_fraction(float stable, float hamers_stable);

float a_out; double M; double f;
double tmax = 2.*M_PI*1.0e2;     // 100 years translated into code units
int n = 0; float inc; float a_over_rh; float mass_ratio = 0.0; float q_in = 1e-2; float stop_inc = 0.0;
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
        a_out = pow(1. + mass_ratio, 1.0/3.0);
        if (argc >= 3) {
            // Only check one given inclination - for computation time
            inc = strtof(argv[2], NULL);
	    if (argc == 4) q_in = strtof(argv[3], NULL);
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

bool hamers_stability_criterion(struct reb_simulation *r){
    r->particles[0].i;
    r->calculate_orbits();
    double final_a_out = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]).a; 
    double final_a_in = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]).a; 
    double a_in = a_over_rh * r_hill * a_out;
    
    bool condition = true;
    if (fabsf(a_in - final_a_in) > 0.1 * fabsf(a_in)) condition = false;
    if (fabsf(a_out - final_a_out) > 0.1 * fabsf(a_out)) condition = false;

    return condition;
}

void loop_sma(){
    /*
     *
    */
    for (a_over_rh = 0.15; a_over_rh < 1.5; a_over_rh += 0.01){
        float stable = 0;
        float hamers_stable = 0;
        for (M = 0; M < 2 * M_PI; M += M_PI/10) {
            f = reb_tools_M_to_f(0, M);
        
            struct reb_simulation* r = reb_create_simulation();
            // Setup constants
            r->dt             = M_PI*1.0e-3;     // initial timestep
            r->integrator        = REB_INTEGRATOR_IAS15;
            r->heartbeat        = heartbeat;

            // Adding star
            struct reb_particle planet = {0}; 
            planet.m  = mass_ratio / (1 + q_in);
            planet.hash = reb_hash("m1");
            reb_add(r, planet); 

            generate_sims(r);

            n++;

            if (r->status != 1) stable++;
            if (hamers_stability_criterion(r)) hamers_stable++;

            reb_free_simulation(r);
        }
        
        print_stable_fraction(stable/20, hamers_stable/20);
    }
}

void generate_sims(struct reb_simulation* r){
    // random mean anomoly
    //f = 2 * M_PI * ((double) random())/ RAND_MAX;
    // the satelite (m2)
    
    struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0],  q_in * r->particles[0].m, a_out * a_over_rh * r_hill, ecc, inc * M_PI / 180, 0.0, 0.0, f);
    m2.hash = reb_hash("m2");
    reb_add(r, m2);

    // The star
    reb_add_fmt(r, "m a e inc", 1.0, a_out, 0.0, 0.0); 

    reb_move_to_com(r);
    
    // add megno particles for stability analysis
    reb_tools_megno_init(r);

    printf("\ninc = %.4f, a = %.4f, n = %d\n", inc, a_over_rh, n);

    char orbits[64];

    if (stop_inc == 0.0){
        sprintf(orbits, "incs/orbits_rand_M%.2e.txt", mass_ratio);
    } else {
        sprintf(orbits, "incs/orbits_rand_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);
    }

    reb_output_orbits(r, orbits);

    reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* r){
    // diagnostics about 1 in every 199 orbits
    double d1 = reb_particle_distance(&r->particles[0], &r->particles[1]);

    if (d1 > a_out/1.5) {
        r->status = 1;
    }
}

void print_megno(struct reb_simulation* r){
    char fname[64];

    sprintf(fname, "incs/megno_rand_M%.2e_i_%.1f.txt", mass_ratio, stop_inc);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.7f\t%d\t%d\t%d\n", reb_tools_calculate_megno(r), n, (int) inc, (int) stop_inc);
    fclose(fp);
}

void print_stable_fraction(float stable, float hamers_stable) {
    char fname[64];
   
    sprintf(fname, "incs/stable_orbits_multiple_M_%.2e_i_%.1f_q_in_%.2e.txt", mass_ratio, stop_inc, q_in);
    
    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.2f\n", stable);
    fclose(fp);

    sprintf(fname, "incs/hamers_stable_orbits_multiple_M_%.2e_i_%.1f_q_in_%.2e.txt", mass_ratio, stop_inc, q_in);
    
    fp = fopen(fname, "a");

    fprintf(fp, "%.2f\n", hamers_stable);
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
