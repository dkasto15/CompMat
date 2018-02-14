// H3 - Albin Annér and David Kasto
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

/* Definitions */
#define nbr_of_timesteps 5000
#define alpha2 0.144
#define h 0.00001
/* Functions */
double V_x(double x);
double weight_function(double timestep, double potential, double energy_trial);
double update_E_T(double energy_trial, double alpha, double timestep,
                  int nbr_of_walkers, int nbr_of_walkers_init);
double V_helium(double x1, double y1, double z1, double x2, double y2, double z2);
double psi_Test(double *r1, double *r2);
double update_walker_1_guided(double * r_walker1, double *r_walker2, int axis, double timestep, double guidance_force, gsl_rng *q);
double update_walker_2_guided(double * r_walker1, double *r_walker2, int axis, double timestep, gsl_rng *q);
double calc_ext_guidance_force_walker_1(double * r_walker1, double *r_walker2, int axis);
double calc_ext_guidance_force_walker_2(double * r_walker1, double *r_walker2, int axis);
double calc_transition_matrix(double x_prev, double x_temp, double timestep, double external_force_prev);

int main() {

	/* Utilility variables */
  FILE *file_energy;
  FILE *file_nbr_of_walkers;
  FILE *file_gauss_test;
  FILE *file_walk_pos1;
  FILE *file_walk_pos2;

  int i, j, k, l, f; /* Looping variables */

	/* Variables associated with random walkers */
	int nbr_of_walkers_init = 300; /* Initial number of walkers */
  int nbr_of_walkers = nbr_of_walkers_init; /* Initial number of walkers */
  double timestep = 0.01; /*0.03;*/ /* Initial time step, should be a value between 0.01
                             and 0.1 */
  double alpha = 0.01; /*0.00001;*/ /*0.002;*/ /* Constant for energy update algorithm, should be a value
                         between 0 and 1 */

  double **pos1 = malloc(nbr_of_walkers*100 * sizeof(double *));
  for(i = 0; i < nbr_of_walkers*100; i++)
    pos1[i] = malloc(3 * sizeof(double)); /* 1D positions for walkers, additional
                                              space left for killing/creating */

  double **pos_temp1 = malloc(nbr_of_walkers*100 * sizeof(double *));
  for(i = 0; i < nbr_of_walkers*100; i++)
    pos_temp1[i] = malloc(3 * sizeof(double)); /* 1D positions for walkers, additional
                                            space left for killing/creating */

  double **pos_prev1 = malloc(nbr_of_walkers*100 * sizeof(double *));
  for(i = 0; i < nbr_of_walkers*100; i++)
    pos_prev1[i] = malloc(3 * sizeof(double)); /* 1D positions for walkers, additional
                                            space left for killing/creating */

  double **pos2 = malloc(nbr_of_walkers *100 * sizeof(double *));
  for(i = 0; i < nbr_of_walkers*100; i++)
    pos2[i] = malloc(3 * sizeof(double)); /* 1D positions for walkers, additional
                                                                                        space left for killing/creating */

  double **pos_temp2 = malloc(nbr_of_walkers*100 * sizeof(double *));
  for(i = 0; i < nbr_of_walkers*100; i++)
    pos_temp2[i] = malloc(3 * sizeof(double)); /* 1D positions for walkers, additional
                                                                                      space left for killing/creating */

  double **pos_prev2 = malloc(nbr_of_walkers*100 * sizeof(double *));
  for(i = 0; i < nbr_of_walkers*100; i++)
      pos_prev2[i] = malloc(3 * sizeof(double)); /* 1D positions for walkers, additional
                                                                                      space left for killing/creating */
  int *nbr_of_walkers_arr = malloc(sizeof(int)*nbr_of_timesteps); /* Array for storing number of walkers
                                                                                                                                                         over simulation procedure */
  double uni_rand; /* Uniform  random number from 0 to 1*/
  double *energy_trial = malloc(sizeof(double)*nbr_of_timesteps); /* Array for storing ground state energy
                                                                                                                                                         oversimulation procedure */
  double gauss_rand; /* Random gaussian number with mean 0 and unit varaince */
  double potential; /* Current potential energy */
  double potential_avg; /* Mean potential energy over time */
  int create_m_walkers; /* Varable for specifying how many walkers that should be created or if the walker
                                                                                                               should be removed */
  int counter; /* Variable for counting how many walkers that have been
                                                                                                      created or removed */

  double guidance_force_prev; /* Guidance force in Force Biased Metropolis */


	/* Initialize uniform random number generator with GSL */
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));
	srand(time(NULL));

	/* * * End of rng initializations * * */

  /* Initialize walkers at random positions */
  for(i = 0; i < nbr_of_walkers_init; i++){
    for(j = 0; j < 3; j++){
      uni_rand = 2*(gsl_rng_uniform(q)-0.5); /* Spread out over -1<x<1 */
      pos_prev1[i][j] = uni_rand;
      uni_rand = 2*(gsl_rng_uniform(q)-0.5); /* Spread out over -1<x<1 */
      pos_prev2[i][j] = uni_rand;
    }
  }
  nbr_of_walkers_arr[0]=nbr_of_walkers_init; /* Initial number of walkers */
  energy_trial[0] = -2.5; /* Initial trial energy */

  file_walk_pos1 = fopen("walk_pos1.dat","w"); /* File for storing the position
                                                  of first array of walkers */
  file_walk_pos2 = fopen("walk_pos2.dat","w"); /* File for storing the position
                                                  of second array of walkers */
    /* SIMULATION */
 	for(k = 1; k < nbr_of_timesteps; k++){
    counter = 0;
    potential_avg=0;

		for(i = 0; i < nbr_of_walkers; i++){

      //for(l = 0; l < 3; l++){
      //  gauss_rand = gsl_ran_gaussian_ratio_method(q,1);
      //  pos_temp1[i][l] = pos_prev1[i][l] + sqrt(timestep)*gauss_rand;
      //  gauss_rand = gsl_ran_gaussian_ratio_method(q,1);
      //  pos_temp2[i][l] = pos_prev2[i][l] + sqrt(timestep)*gauss_rand;
      //}

      /* Update particles according to Force Biased Diffusion Monte Carlo method */
      for(l = 0; l < 3; l++){
        guidance_force_prev = calc_ext_guidance_force_walker_1(pos_prev1[i], pos_prev2[i], l);

        pos_temp1[i][l] = update_walker_1_guided(pos_prev1[i], pos_prev2[i], l, timestep, guidance_force_prev, q);

        guidance_force_prev = calc_ext_guidance_force_walker_2(pos_prev1[i], pos_prev2[i], l);

        pos_temp2[i][l] = update_walker_2_guided(pos_prev1[i], pos_prev2[i], l, timestep, q);
      }

    	potential = V_helium(pos_temp1[i][0], pos_temp1[i][1], pos_temp1[i][2],
                           pos_temp2[i][0], pos_temp2[i][1], pos_temp2[i][2]);
      potential_avg += V_helium(pos_prev1[i][0], pos_prev1[i][1], pos_prev1[i][2],
                           pos_prev2[i][0], pos_prev2[i][1], pos_prev2[i][2]);
    	uni_rand = gsl_rng_uniform(q);

			create_m_walkers = (int) (weight_function(timestep, potential,
                         energy_trial[k-1]) + uni_rand);

			if(create_m_walkers > 0){
				for(j = 0; j < create_m_walkers; j++){
          for(l = 0; l < 3; l++){
					  pos1[counter][l] = pos_temp1[i][l];
            pos2[counter][l] = pos_temp2[i][l];
          }
          if (k > (nbr_of_timesteps - 20000)) {
            fprintf(file_walk_pos1, "%e %c %e %c %e %c \n", pos_temp1[i][0], ',',
                                                            pos_temp1[i][1], ',',
                                                            pos_temp1[i][2], ',');
            fprintf(file_walk_pos2, "%e %c %e %c %e %c \n", pos_temp2[i][0], ',',
                                                            pos_temp2[i][1], ',',
                                                            pos_temp2[i][2], ',');
          }
          counter += 1;
				}
			} else if(create_m_walkers == 0) {
        ;
			} else {
				printf("%s\n", "Something went wrong: Creating negative number of walkers");
			}
    }

    potential_avg = potential_avg/((double)nbr_of_walkers);
    nbr_of_walkers = counter;

  	energy_trial[k] = update_E_T(potential_avg, alpha, timestep,
                                 nbr_of_walkers, nbr_of_walkers_init);

    if(fabs(energy_trial[k]) > 7){
      printf("%s\n","**** LOW ENERGY ****" );
      printf("%s%e\n", "E_T: \t", energy_trial[k]);
      printf("%s%e\n","Average potential energy: \t", potential_avg);
      printf("%s%e\n","Log number of walkers quota: \t", - (alpha/timestep)
             *log((double)nbr_of_walkers/(double)nbr_of_walkers_init));
      printf("%s\n", "********");
    }

    nbr_of_walkers_arr[k] = nbr_of_walkers;

    for (i = 0; i < nbr_of_walkers; i++) {
      for(l = 0; l < 3; l++){
        pos_prev1[i][l]=pos1[i][l];
        pos_prev2[i][l]=pos2[i][l];
      }
    }
    printf("%d\n", nbr_of_walkers);
	}

  file_energy = fopen("energy.dat", "w");
  file_nbr_of_walkers = fopen("nbr_of_walkers.dat", "w");

  for(k = 0; k < nbr_of_timesteps; k++){
    fprintf(file_energy, "%e \n", energy_trial[k]);
    fprintf(file_nbr_of_walkers, "%d \n", nbr_of_walkers_arr[k]);
  }

  fclose(file_energy);
  fclose(file_nbr_of_walkers);
  fclose(file_walk_pos1);
  fclose(file_walk_pos2);

  for(i = 0; i < nbr_of_walkers*100; i++)
    free((void*)pos1[i])
    free((void*)pos_temp1[i])
    free((void*)pos_prev1[i])
    free((void*)pos2[i])
    free((void*)pos_temp2[i])
    free((void*)pos_prev2[i])

  free(nbr_of_walkers_arr)
  free(energy_trial)

}

double V_x(double x){
  return 0.5*x*x;
}

double V_helium(double x1, double y1, double z1, double x2, double y2, double z2){
  return -2.0/sqrt(x1*x1 + y1*y1 + z1*z1) -2.0/sqrt(x2*x2 + y2*y2 + z2*z2)
          + 1/sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

double weight_function(double timestep, double potential, double energy_trial){
  return exp(-timestep*(potential-energy_trial));
}

double update_E_T(double energy_trial, double alpha, double timestep,
                  int nbr_of_walkers, int nbr_of_walkers_init){
  return energy_trial - (alpha/timestep)
         *log((double)nbr_of_walkers/(double)nbr_of_walkers_init);
}
double psi_Test(double *r1, double *r2){
	double r1_abs = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
	double r2_abs = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
	double r12_abs = sqrt((r1[0]-r2[0])*(r1[0]-r2[0])+
										(r1[1]-r2[1])*(r1[1]-r2[1])+
										(r1[2]-r2[2])*(r1[2]-r2[2]));
	double psi_T = exp(-2.0*r1_abs)*exp(-2.0*r2_abs)*exp(r12_abs/(2.0*(1.0+alpha2*r12_abs)));
	return psi_T;
}

double update_walker_1_guided(double * r_walker1, double *r_walker2, int axis, double timestep, double guidance_force, gsl_rng *q){
  return r_walker1[axis] +
         timestep/2*guidance_force +
         sqrt(timestep)*gsl_ran_gaussian_ratio_method(q,1);
}

double update_walker_2_guided(double * r_walker1, double *r_walker2, int axis, double timestep, gsl_rng *q){
  double guidance_force;

  guidance_force = calc_ext_guidance_force_walker_2(r_walker1, r_walker2, axis);
  return r_walker1[axis] +
         timestep/2*guidance_force +
         sqrt(timestep)*gsl_ran_gaussian_ratio_method(q,1);
};

double calc_ext_guidance_force_walker_1(double * r_walker1, double *r_walker2, int axis){
  double F;
  double dpsi_dr_walker;
  double r_plus_h[3];
  double r_minus_h[3];
  double psi_plus_h;
  double psi_minus_h;
  double psi = psi_Test(r_walker1, r_walker2);

  int i;
  for(i = 0; i < 3; i++){
    if(i == axis){
      r_plus_h[i] = r_walker1[i] + h;
      r_minus_h[i] = r_walker1[i] - h;
    } else {
      r_plus_h[i] = r_walker1[i];
      r_minus_h[i] = r_walker1[i];
    }
  }

  psi_plus_h = psi_Test(r_plus_h, r_walker2);
  psi_minus_h = psi_Test(r_minus_h, r_walker2);
  dpsi_dr_walker = (psi_plus_h - 2.0*psi + psi_minus_h)/(h*h);

  F = (2.0/psi)*dpsi_dr_walker;
  return F;
}

double calc_ext_guidance_force_walker_2(double * r_walker1, double *r_walker2, int axis){
  double F;
  double dpsi_dr_walker;
  double r_plus_h[3];
  double r_minus_h[3];
  double psi_plus_h;
  double psi_minus_h;
  double psi = psi_Test(r_walker1, r_walker2);

  int i;
  for(i = 0; i < 3; i++){
    if(i == axis){
      r_plus_h[i] = r_walker2[i] + h;
      r_minus_h[i] = r_walker2[i] - h;
    } else {
      r_plus_h[i] = r_walker2[i];
      r_minus_h[i] = r_walker2[i];
    }
  }

  psi_plus_h = psi_Test(r_walker1, r_plus_h);
  psi_minus_h = psi_Test(r_walker1, r_minus_h);
  dpsi_dr_walker = (psi_plus_h - 2.0*psi + psi_minus_h)/(h*h);

  F = (2.0/psi)*dpsi_dr_walker;
  return F;
}
