#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/resource.h>
#include "Gradient_Sensing_Cell_ML_git_stamps.h"



/* This code generates N trajectories of Brownian particles.
 * They eminate from a source at the origin and proceeed until
 * their distances is greater than cutoff or until they hit the
 * surface of a sphere (center on z-axis at distance d radius r).
 *
 * The output is 
 * time coordinate
 *
 * Plus other output. 
 *
 * Potential optimisation: Make every trajectory count by assuming that 
 * a sphere was in the way of a trajectory that reaches the cutoff.
 *
 * At the moment, I don't bother much about optimising -- I draw normally
 * distributed rvs.
 *
 *
 * make BrownianParticle
 * ./BrownianParticle > BrownianParticle_ref.dat &
 * grep EVENT BrownianParticle_ref.dat | sed 's/.*EVENT //' > BrownianParticle_ref.txt
 * gnuplot> sp 'BrownianParticle_ref.txt' u 5:6:7
 *
 * ./BrownianParticle -p 1.001 > BrownianParticle_ref2.dat
 * ... has a large number of particles arriving at the point closest to the origin.
 *
 * Validation:
 * Look at the histogram of column 5, that's the x-coordinate at arrival.
 * How does that compare to theory?
 */


/*
 * double gsl_ran_gaussian_ziggurat(const gsl_rng *r, double sigma)
 */


/* Distance the cue particle has to diffusive before considered "lost". */
#ifndef CUTOFF
#define CUTOFF (20.)
#endif
double param_cutoff=CUTOFF;
double param_cutoff_squared=CUTOFF*CUTOFF;

/* Radius of the cell sphere. */ 
#ifndef SPHERE_RADIUS
#define SPHERE_RADIUS (1.0)
#endif
double param_sphere_radius=SPHERE_RADIUS;
double param_sphere_radius_squared=SPHERE_RADIUS*SPHERE_RADIUS;

/* Position of the cell sphere along x-axis. */
#ifndef SPHERE_POSZ
#define SPHERE_POSZ (5.0)
#endif
double param_sphere_posz=SPHERE_POSZ;


/* Diffusive step length, sigma=2D\Delta t in every direction. */
#ifndef SIGMA 
#define SIGMA (0.01)
#endif
double param_sigma=SIGMA;

/* Number of iterations */
#ifndef NUM_ITERATIONS
#define NUM_ITERATIONS (1000000LL)
#endif
long long int param_num_iterations=NUM_ITERATIONS;
long long int param_num_events=0LL;

/* Seed */
#ifndef SEED
#define SEED (5UL)
#endif
unsigned long int param_seed=SEED;

gsl_rng *rng;



void postamble(void);


int main(int argc, char *argv[])
{
double x, y, z;
long long int iteration;
long long int event=0LL;
int ch;
double origin_distance2; /* Variable names with a 2 at the end are squared, parameters have that suffix, to avoid problems by typos. */
double sphere_distance2, distance_from_z_axis2;
long long int steps;

setlinebuf(stdout);

{ 
int i;

printf("# Info: Command: %s", argv[0]);
for (i=1; i<argc; i++) printf(" \"%s\"", argv[i]);
printf("\n");
}


/* Some infos. */
{
time_t tim;
tim=time(NULL);

printf("# Info Starting at %s", ctime(&tim));
}

/* For version control if present. */
printf("# Info: Version of git_version_string to follow.\n");
printf("%s", git_version_string);
printf("# $Header$\n");


/* Hostname */
{ char hostname[128];
  gethostname(hostname, sizeof(hostname)-1);
  hostname[sizeof(hostname)-1]=(char)0;
  printf("# Info: Hostname: %s\n", hostname);
}


/* Dirname */
{ char cwd[1024];
  cwd[0]=0;
  if(getcwd(cwd, sizeof(cwd)-1)!=NULL){
  cwd[sizeof(cwd)-1]=(char)0;
  printf("# Info: Directory: %s\n", cwd);}
}

/* Process ID. */
printf("# Info: PID: %i\n", (int)getpid());

while ((ch = getopt(argc, argv, "c:E:N:p:r:s:S:")) != -1) {
  switch (ch) {
    case'c':
      param_cutoff=strtod(optarg, NULL);
      break;
    case 'E':
      param_num_events=strtoll(optarg, NULL, 10);
      break;
    case 'N':
      param_num_iterations=strtoll(optarg, NULL, 10);
      break;
    case 'p':
      param_sphere_posz=strtod(optarg, NULL);
      break;
    case 's':
      param_sigma=strtod(optarg, NULL);
      break;
    case 'S':
      param_seed=strtoul(optarg, NULL, 10);
      break;
    default:
      break;
    }
  }

rng=gsl_rng_alloc(gsl_rng_taus2);
gsl_rng_set(rng, SEED);


#define PRINT_PARAM(a,o, f) printf("#Info: %s: %s " f "\n", #a, o, a)


param_cutoff_squared=param_cutoff*param_cutoff;
param_sphere_radius_squared=param_sphere_radius*param_sphere_radius;
PRINT_PARAM(param_num_events, "-E", "%lli");
PRINT_PARAM(param_num_iterations, "-N", "%lli");
PRINT_PARAM(param_sigma, "-s", "%g");
PRINT_PARAM(param_cutoff, "-c", "%g");
PRINT_PARAM(param_cutoff_squared, "", "%g");
PRINT_PARAM(param_sphere_radius, "", "%g");
PRINT_PARAM(param_sphere_radius_squared, "", "%g");
PRINT_PARAM(param_sphere_posz, "-p", "%g");
PRINT_PARAM(param_seed, "-S", "%lu");



for (iteration=1LL; ((iteration<=param_num_iterations) &&   ( (param_num_events==0) ? (1) : (event<param_num_events))); iteration++) {
  x=y=z=0.;
  origin_distance2=0.;
  steps=0LL;

  while (origin_distance2<=param_cutoff_squared) {
    steps++;
    x+=gsl_ran_gaussian_ziggurat(rng, param_sigma);
    y+=gsl_ran_gaussian_ziggurat(rng, param_sigma);
    z+=gsl_ran_gaussian_ziggurat(rng, param_sigma);



    origin_distance2 = sphere_distance2 = distance_from_z_axis2 = y*y+x*x;
    origin_distance2+=z*z;

    sphere_distance2+=(z-param_sphere_posz)*(z-param_sphere_posz);

    //printf("Particle %lli event %lli step %lli pos %g %g %g distances %g %g\n", iteration, event, steps, x, y, z, origin_distance2, sphere_distance2);
    if (sphere_distance2<param_sphere_radius_squared) {
      double theta, phi;

      /* Shift x by param_sphere_posx and then convert coordinates to 
       * polar and azimuthal angles. 
       * The azimuthal angle goes 0...\pi,
       * but atan returns 0..\pi/2 for 0..\infty
       * and then -\pi/2..0 for -\infty..-0 */
      z=param_sphere_posz-z; /* Distance from sphere centre. */
      if (z!=0.) {
        if ((theta=atan(sqrt(distance_from_z_axis2)/z))<0.) theta+=M_PI;
      } else theta=M_PI/2.;
      phi=atan2(y,x); /* phi=0 for y=0 */

      event++;

      printf("# EVENT %lli %lli %10.20g %10.20g %10.20g %10.20g %10.20g %lli\n", event, iteration, theta, phi, x, y, z, steps);
      /* Watch out, z is now relative to the origin of the sphere. It has been shifted by param_sphere_posz above. */
      break;
    }
  }
  if (origin_distance2>param_cutoff_squared)
    printf("# PARTICLE %lli got lost to position %g %g %g after %lli steps at distance %g>%g.\n", iteration, x, y, z, steps, sqrt(origin_distance2), param_cutoff);
}

postamble();
return(0);
}






void postamble(void)
{
time_t tm;
struct rusage rus;

tm=time(NULL);

printf("# Info Terminating at %s", ctime(&tm));
if (getrusage(RUSAGE_SELF, &rus)) {
  printf("# Info getrusage(2) failed.\n");
} else {
printf("# Info getrusage.ru_utime: %li.%06li\n", (long int)rus.ru_utime.tv_sec, (long int)rus.ru_utime.tv_usec);
printf("# Info getrusage.ru_stime: %li.%06li\n", (long int)rus.ru_stime.tv_sec, (long int)rus.ru_stime.tv_usec);

#define GETRUSAGE_long(f) printf("# Info getrusage.%s: %li\n", #f, rus.f);
GETRUSAGE_long(ru_maxrss);
GETRUSAGE_long(ru_ixrss);
GETRUSAGE_long(ru_idrss);
GETRUSAGE_long(ru_isrss);
GETRUSAGE_long(ru_minflt);
GETRUSAGE_long(ru_majflt);
GETRUSAGE_long(ru_nswap);
GETRUSAGE_long(ru_inblock);
GETRUSAGE_long(ru_oublock);
GETRUSAGE_long(ru_msgsnd);
GETRUSAGE_long(ru_msgrcv);
GETRUSAGE_long(ru_nsignals);
GETRUSAGE_long(ru_nvcsw);
GETRUSAGE_long(ru_nivcsw);
}
printf("#Info: Good bye and thanks for all the fish.\n");
}

