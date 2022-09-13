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
void make_orbits();
bool stopped;
double tmax = 2.*M_PI*1.0e1;     // 100 years translated into code units
int n = 0; float inc; float a_over_rh; float q = 0.0; float q_in = 0;
float ecc = 0.0, a_out;
// For plotting info about a specific orbit
int target_inc = 50;
float target_a = 0.4, e_in, e_out, M1, M2, f1, f2, omeg1, Omeg1, a_over_rp, a_in;
long int seed;

// hill radius
float r_hill;

int main(int argc, char* argv[]){
    // Initial conditions
    if (argc == 1) seed = 1;
    else seed = strtol(argv[1], NULL, 10);
    srandom(seed);
    make_orbits();
}

float predict_stable(float inc, float q){
    // Predict if an orbit is stable given its orbital parameters
    float a, b, stable_a_over_rp;
    if (log10(q) <= -0.41){
        a = 0.32;
        b = -0.52;
    } else {
        a = 0.42;
        b = -0.48;
    }
    float q_factor = pow(10, b) * pow(q, a);
    float combination_factor = pow(10, (inc * pow(q, 1.3))/ 1500);
    
    double inc_rad = inc * M_PI / 180;
    float inc_factor;
    if (inc < 60) {
        inc_factor = -0.17 * cos(inc_rad) + 0.63;
    }
    else {
        float a0 = -0.0798, a1 = 0.5045, a2 = -0.8617, a3 = 0.2260, a4 = 0.7300;
        inc_factor = a0 * pow(inc_rad, 4) + a1 * pow(inc_rad, 3) + a2 * pow(inc_rad, 2) + a3 * inc_rad + a4;
    }
    inc_factor /= 0.45;
	
    printf("inc factor = %.5f", inc_factor);
    printf("q factor = %.5f", q_factor);
    printf("combination factor = %.5f", combination_factor);


    stable_a_over_rp = q_factor * inc_factor / combination_factor;
    
    return stable_a_over_rp; 
}

double mardling(double q, double inc, double e_out){
    double non_inc = pow(((1 + 1 / q) * (1 + e_out)/(sqrt(1 - e_out))), -0.4) / 2.8;
    return non_inc / (1 - 0.3 * inc / 180);
}

void make_orbits(){
    /*
     *
    */
        
    for (int i = 0; i < 4000; i++){
        inc = 180 * ((double) random()) / RAND_MAX;
	    double a = pow(10, -2 * ((double) random()) / RAND_MAX), 
		b = pow(10, -2 * ((double) random()) / RAND_MAX), 
		c = pow(10, -2 * ((double) random()) / RAND_MAX), 
        q = a / (b + c);
    	q_in = b / c < c / b ? b/c : c/b;
        r_hill = pow(q/3, 1./3.);
        e_in = ((double) random()) / RAND_MAX;
        e_out = ((double) random()) / RAND_MAX;
        M1 = 2 * M_PI * ((double) random()) / RAND_MAX;
        M2 = 2 * M_PI * ((double) random()) / RAND_MAX;
        Omeg1 = 2 * M_PI * ((double) random()) / RAND_MAX;
        omeg1 = 2 * M_PI * ((double) random()) / RAND_MAX;
        a_over_rp = (0.9999 * ((double) random()) / RAND_MAX + 1e-4) / (1 - e_out);
        a_out = pow(1 + q, 1./3.);

        struct reb_simulation* r = reb_create_simulation();
        // Setup constants
        r->dt             = M_PI*1.0e-3;     // initial timestep
        r->integrator        = REB_INTEGRATOR_IAS15;
        r->heartbeat        = heartbeat;
            // Adding star
        struct reb_particle planet = {0}; 
        planet.m  = q/(1 + q_in);
        planet.hash = reb_hash("m1");
        printf("m1 = %d", planet.hash);
        reb_add(r, planet); 

        generate_sims(r);
        
        // checking for >10% change in a_in at the end of the 100 periods
        char fname[64];

        sprintf(fname, "fit_accuracy/stability_predictions_seed_%d.txt", seed);

        FILE* fp = fopen(fname, "a");
        fprintf(fp, "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\n", log10(q), log10(q_in), inc, e_in, e_out, a_over_rp, r->status!=1);
        fclose(fp);

        n++;

        reb_free_simulation(r);
    }
}

void generate_sims(struct reb_simulation* r){
    // translating mean anomalies to true anomalies
    f1 = reb_tools_M_to_f(e_in, M1);
    f2 = reb_tools_M_to_f(e_out, M2);

    // the smaller inner body
    struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0], q_in * r->particles[0].m, a_over_rp * (1 - e_out) * a_out, e_in, inc * M_PI / 180, Omeg1, omeg1, f1);
    m2.hash = reb_hash("m2");
    printf("\nm2 = %d", m2.hash);
    reb_add(r, m2);

    // The star
    struct reb_particle star = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0, a_out, e_out, 0.0, 0.0, 0.0, f2);
    star.hash = reb_hash("star"); 
    printf("\nstar = %d", star.hash);
    reb_add(r, star); 

    reb_move_to_com(r);
    
    printf("\ninc = %.4f, a = %.4f, n = %d\n", inc, a_over_rp, n);
    stopped = false;
    reb_integrate(r, tmax);
    // when orbit was either finished or disrupted, output time
}

void heartbeat(struct reb_simulation* r){
	float a_in = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]).a;
	float current_a_out = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]).a; // difficult to calculate accurate a_out when q_in ~ 1
	float initial_a_in = a_over_rp * (1 - e_out) * a_out;
	// diagnostics about 1 in every 199 orbits

	// exiting when a_in or a_out changes by more than 10%
	if (fabsf(a_in - initial_a_in) > 0.1 * fabsf(initial_a_in) || fabsf(current_a_out - a_out) > 0.1 * fabsf(a_out)){
		r->status = 1;
	}
}

