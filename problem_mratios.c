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
void print_stable_fraction(float stable);

double tmax = 2.*M_PI*1.0e2;     // 100 years translated into code units
int n = 0; float inc; float a_over_rh; float mass_ratio = 0.0; double f, M;
float ecc = 0.0, mass_ratio_exp = 0.0, stop_mass_ratio_exp = 0.0;
// hill radius
float a_out, ain_over_aout;

int main(int argc, char* argv[]){
    // Initial conditions
    if (argc > 1){
        // Mass ratio between planet and outer star
        mass_ratio_exp = strtof(argv[1], NULL);


        if (argc == 4) {
            // Only check one given inclination - for computation time
            inc = strtof(argv[2], NULL);
            ecc = strtof(argv[3], NULL);
            stop_mass_ratio_exp = mass_ratio_exp - 0.1;    
            //srandom((int) (mass_ratio_exp * 1000));    
            for (; mass_ratio_exp > stop_mass_ratio_exp; mass_ratio_exp -=0.01) {
		mass_ratio = pow(10, mass_ratio_exp);
		// to preserve period of 2*pi time units (kepler 3rd law)
        	a_out = pow((mass_ratio + 1.0), 1.0/3.0);
                loop_sma();
            }
        } 
        else {
            printf("wrong number of command line arguments");
        }
    }
    else {
        printf("Must provide mass ratio as command line argument!\n");
        return 1;
    }    
}

bool approx_stable_criterion(double log_ain_over_aout, double log_q){
    return (log_ain_over_aout < (log_q/3 - 0.65));
}

bool approx_unstable_criterion(double log_ain_over_aout, double log_q){
    return (log_ain_over_aout > (log_q/3 - 0.1));
}

void loop_sma(){
    /*
     *
    */
    for (double log_ain_over_aout = -1.80; log_ain_over_aout <= -0.299; log_ain_over_aout += 0.01){
        float stable = 0;	
        
        if (approx_stable_criterion(log_ain_over_aout, mass_ratio_exp)) {
            print_stable_fraction(1);
            continue;
        }

        if (approx_unstable_criterion(log_ain_over_aout, mass_ratio_exp)) {
            print_stable_fraction(0);
            continue;
        }

        for (M = 0; M < 2 * M_PI; M += M_PI / 10) {
	    f = reb_tools_M_to_f(0, M);

            ain_over_aout = pow(10, log_ain_over_aout);
            
            struct reb_simulation* r = reb_create_simulation();
            // Setup constants
            r->dt             = M_PI*1.0e-6;     // initial timestep
            r->integrator        = REB_INTEGRATOR_IAS15;
            r->heartbeat        = heartbeat;

            // Adding star
            struct reb_particle planet = {0}; 
            planet.m  = mass_ratio;
            planet.hash = reb_hash("m1");
            reb_add(r, planet); 

            generate_sims(r);

            n++;

            if (r->status != 1) stable++;

            reb_free_simulation(r);
        }
        
        print_stable_fraction(stable/20);       
    }
}

void generate_sims(struct reb_simulation* r){
    // pre-determined mean anomaly
    //f = (((double) random())/RAND_MAX) * 2 * M_PI; 
    struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0e-4 * r->particles[0].m, ain_over_aout * a_out, ecc, inc * M_PI / 180, 0.0, 0.0, f);
    m2.hash = reb_hash("m2");
    reb_add(r, m2);

    // The star
    struct reb_particle star = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0, a_out, 0.0, 0.0, 0.0, 0.0, 0.0);
    star.hash = reb_hash("star"); 

    reb_add(r, star); 

    reb_move_to_com(r);
    
    printf("\na = %.4f, q = %.6f, n = %d\n", ain_over_aout, mass_ratio, n);

    char orbits[64];

    sprintf(orbits, "mass_ratios/orbits_0_M_%.2f_i_%.1f.txt", stop_mass_ratio_exp, inc);

    reb_output_orbits(r, orbits);

    reb_integrate(r, tmax);
    // when orbit was either finished or disrupted, output time
    print_time(r);
}

void heartbeat(struct reb_simulation* r){
    // distance from the satellite to the earth vs the sun
    float d_in = reb_particle_distance(&r->particles[1], &r->particles[0]);
    double e_in = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]).e;
    if(d_in > a_out/1.5 || e_in > 0.99){          
        r->status = 1;
    }
}

void print_time(struct reb_simulation* r){
    char fname[64];

    sprintf(fname, "mass_ratios/time_0_M_%.2f_i_%.2f.txt", stop_mass_ratio_exp, inc);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.7f\t%d\t%.2f\t%.2f\n", r->t, n, mass_ratio_exp, stop_mass_ratio_exp);
    fclose(fp);
}

void print_stable_fraction(float stable) {
    char fname[64];

    sprintf(fname, "mass_ratios/stable_orbits_multiple_M_%.2f_i_%.1f.txt", stop_mass_ratio_exp, inc);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.2f\n", stable);
    fclose(fp);
}

