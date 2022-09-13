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
void print_stable(float stable_frac);
void generate_sims(struct reb_simulation* r);
void global_orbit_output(struct reb_simulation* r, struct reb_particle p);
void relative_orbit_output(struct reb_simulation* r, struct reb_particle s, struct reb_particle p);
void loop_sma();

double tmax = 2.*M_PI*1.0e2;     // 100 years translated into code units
int n = 0; float inc; float a_over_rhp; float mass_ratio = 0.0; double M;
float ecc = 0.0, stop_ecc = 0.0;
// For plotting info about a specific orbit
int target_inc = 50;
float target_a = 0.4;
double a_out;

double meg = 0;

// hill radius
float r_hill;

int main(int argc, char* argv[]){
    // Initial conditions
    printf("Started main function!\n");
    if (argc > 1){
        // Mass ratio between planet and outer star
        mass_ratio = pow(10, strtof(argv[1], NULL));

        r_hill = pow((mass_ratio)/3., 1.0/3.0);
        a_out = pow(mass_ratio + 1, 1.0/3.0);

        if (argc == 4) {
            // Only check one given inclination - for computation time
            //inc = strtof(argv[2], NULL);
            // Random inc
            inc = strtof(argv[2], NULL);
            inc *= 180/M_PI;
            ecc = strtof(argv[3], NULL);
            stop_ecc = ecc + 0.05;    
            printf("mrat, inc, ecc = %.3f, %.3f, %.3f", mass_ratio, inc, ecc);
            while (ecc < stop_ecc) {
                loop_sma();
                ecc += 0.01;
            }
        } 
	else if (argc == 3){
	    printf("wrong number of command line arguments: received %d", argc);
	    printf("args were %s %s %s", argv[0], argv[1], argv[2]);
	}
    }
    else {
        printf("Must provide mass ratio as command line argument!\n");
        return 1;
    }    
}

float predict_stable(float inc, float mass_ratio){
    // Predict if an orbit is stable given its orbital parameters
    float a, b, stable_a_over_rp;
    if (log10(mass_ratio) <= -0.41){
        a = 0.32;
        b = -0.52;
    } else {
        a = 0.42;
        b = -0.48;
    }

    float mass_ratio_factor = pow(10, b) * pow(mass_ratio, a);
    float combination_factor = pow(10, (inc * pow(mass_ratio, 1.3))/ 1500);

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

    stable_a_over_rp = mass_ratio_factor * inc_factor / combination_factor;

    return stable_a_over_rp;
}


void loop_sma(){
    /*
     *
    */
    for (a_over_rhp = 0.2; a_over_rhp <= 1.0; a_over_rhp += 0.01){
        float stable = 0;
        static int unstable = 0;
        if (a_over_rhp <= 0.3) {
            unstable = 0;
            printf("resetting to 0 unstable\n");
        }

        if (unstable >= 5){
            printf("all unstable");
            print_stable(0.0);
            continue;
        }

        // eliminating clear cases to save time
        if (a_over_rhp < predict_stable(inc, mass_ratio) / r_hill - 0.1){
            printf("all stable");
            print_stable(1.0);
            continue;
        }
        if (a_over_rhp > predict_stable(inc, mass_ratio) / r_hill + 0.3){
            printf("all unstable");
            print_stable(0.0);
            continue;
        }
        
        for (M = 0; M < 2*M_PI; M += M_PI/10){

            struct reb_simulation* r = reb_create_simulation();
            // Setup constants
            r->dt             = M_PI*1.0e-6;     // initial timestep
            r->integrator        = REB_INTEGRATOR_IAS15;
            r->heartbeat        = heartbeat;
            r->collision = REB_COLLISION_DIRECT;

            // Adding planet
            struct reb_particle planet = {0}; 
            planet.m  = mass_ratio;
            planet.hash = reb_hash("m1");
            reb_add(r, planet); 

            generate_sims(r);

            if (r->status != 1)stable++;
            n++;

            reb_free_simulation(r);
        }
        
        print_stable(stable/20);

        if (stable == 0){
            unstable++;
        }
    }
}

void generate_sims(struct reb_simulation* r){
    double f = reb_tools_M_to_f(ecc, M);
    
    // The moon (m2)
    double r_hp = r_hill * (1 - ecc);
    struct reb_particle m2 = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0e-2 * r->particles[0].m, a_over_rhp * r_hp, 0, inc * M_PI / 180, 0.0, 0.0, f);
    m2.hash = reb_hash("m2");
    reb_add(r, m2);

    // The star
    // Initialise at true anomaly of pi - this puts the star at apoapsis initially
    struct reb_particle star = reb_tools_orbit_to_particle(r->G, r->particles[0],  1.0, a_out, ecc, 0.0, 0.0, 0.0, M_PI);
    star.hash = reb_hash("star"); 

    reb_add(r, star); 

    reb_move_to_com(r);

    printf("\ninc = %.4f, a = %.4f, n = %d\n", inc, a_over_rhp, n);

    char orbits[64];

    sprintf(orbits, "eccentricities/orbits_0_M_%.2e.txt", mass_ratio);

    reb_output_orbits(r, orbits);

    reb_integrate(r, tmax);
    // when orbit was either finished or disrupted, output time
    print_time(r);
}

void heartbeat(struct reb_simulation* r){
    // diagnostics about 1 in every 199 orbits
    double d_in = reb_particle_distance(&r->particles[0], &r->particles[1]);
    double e = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]).e;

    if (d_in > a_out / 1.5 || e > 0.99) {
        r->status = 1;
    }
}

void print_stable(float stable_frac){
    char fname[64];

    sprintf(fname, "eccentricities_out/stable_orbits_%.2e_e_%.2f.txt", mass_ratio, stop_ecc);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.2f\n", stable_frac);
    fclose(fp);
}

void print_time(struct reb_simulation* r){
    char fname[64];

    sprintf(fname, "eccentricities_out/time_0_M_%.2e_e_%.2f.txt", mass_ratio, stop_ecc);

    FILE* fp = fopen(fname, "a");

    fprintf(fp, "%.7f\t%d\t%.2f\t%.2f\n", r->t, n, ecc, stop_ecc);
    fclose(fp);
}
