// for geneerating 
//derived from the REBound example
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

#include "hdf5.h"
#include "hdf5_hl.h"

// disc model
const double Sigma0       = 3.8e-4 *0.5*1.4;
const double redge        = 0.1;
const double deltaredge   = 0.001;
const double aspectratio0 = 0.035;
const double alpha        = 1.0;
const double flaringindex = 2.0/7.0;
const double ffudge       = 1.0/100.0; //fudge factor for matsumoto et al. 2012 - not sure on sign

const double tmax      = 16.0e4*2.*M_PI;    // in year/(2*pi)
const double tdep      = 0.3*tmax;
const double deltatdep = 0.2*tmax;
const double dtout     = 40.0;

const double starmass = 1.0;

//parameters for IC
const int p = 4;
const int q = 1;
const int nchain = 16; 
const double pmass = 1e-5;
double ai = 0.11;

const int noutmax = (int)(tmax/dtout) +1;
int iout;

struct output_structure {
  int nout;
  int nchain;
  int p;
  int q;
  double G;
  double starmass;
  double *pmass; 
  double ai;
  double Sigma0;
  double redge;
  double deltaredge;
  double aspectratio0;
  double alpha;
  double flaringindex;
  double ffudge;
  double tmax;
  double tdep;
  double deltatdep;
  double *t;
  // the memebrs of the reb_orbit structure, as arrays
  double *d;        ///< Radial distance from central object
  double *v;        ///< velocity relative to central object's velocity
  double *h;        ///< Angular momentum
  double *P;        ///< Orbital period
  double *n;        ///< Mean motion
  double *a;        ///< Semi-major axis
  double *e;        ///< Eccentricity
  double *inc;      ///< Inclination
  double *Omega;    ///< Longitude of ascending node
  double *omega;    ///< Argument of pericenter
  double *pomega;   ///< Longitude of pericenter
  double *f;        ///< True anomaly
  double *M;        ///< Mean anomaly
  double *l;        ///< Mean Longitude
  double *theta;    ///< True Longitude
  double *T;        ///< Time of pericenter passage
  double *rhill;    ///< Circular Hill radius 
  double *Sigma;    ///< Disc gas surface denisty at planet location
};

struct output_structure out;

void migration_forces(struct reb_simulation* r);
void migration_forces_tanaka_ward(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);
void init_output_structure(struct output_structure* out, int nchain, int nout);
void out_to_hdf5(struct output_structure* out);
void free_output_structure(struct output_structure* out);
double disc_surface_density(const double r, const double t);
double get_particle_r(struct reb_particle* p, struct reb_particle* com);

int main(int argc, char* argv[]){
    // set output counter
    iout = 0;
    init_output_structure(&out, nchain, noutmax);

    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->integrator    = REB_INTEGRATOR_WHFAST;
    //r->integrator    = REB_INTEGRATOR_IAS15;
    r->dt         = 1e-2*2.*M_PI *pow(0.1, 1.5);        // in year/(2*pi)
    r->additional_forces = migration_forces_tanaka_ward;     //Set function pointer to add dissipative forces.
    r->heartbeat = heartbeat;  
    r->force_is_velocity_dependent = 1;

    // Initial conditions
    // Parameters are those of Lee & Peale 2002, Figure 4. 
    struct reb_particle star = {0};
    star.m  = out.starmass;            // This is a sub-solar mass star
    reb_add(r, star); 
    
    int errcode;
    for (int ip=0; ip < nchain; ip++){
      // struct reb_particle reb_tools_orbit_to_particle_err(double G, struct reb_particle primary, 
      //      double m, double a, double e, double i, double Omega, double omega, double f, int* err)
      struct reb_particle pt = reb_tools_orbit_to_particle_err( 
              r->G, star, pmass, ai, ip*0.01, ip*0.01, 0.0, 0.0, 0.0, &errcode);
      reb_add(r, pt);
      ai = ai * pow(p / (p + q + 0.05), -2./3.);
    } 


    reb_move_to_com(r);          

    system("rm -v orbits.txt"); // delete previous output file

    reb_integrate(r, tmax);
   
    out_to_hdf5(&out);
    free_output_structure(&out);

    printf("\n");
}

void init_output_structure(struct output_structure* out, int nchain, int nout){
    out->nout = nout;
    out->nchain = nchain;
    out->tmax = tmax;
    out->tdep = tdep;
    out->deltatdep = deltatdep;
    out->Sigma0 = Sigma0;
    out->redge  = redge;
    out->deltaredge = deltaredge;
    out->aspectratio0 = aspectratio0;
    out->alpha = alpha;
    out->flaringindex = flaringindex;
    out->ffudge = ffudge;
    out->G = 1.0;
    out->starmass = starmass;
  
    out->p = p;
    out->q = q;
    out->pmass = calloc(nchain, sizeof(double));;
    for (int i=0; i<nchain; i++)
        out->pmass[i] = pmass;
    out->ai = ai;

    out->t = calloc(nchain*nout, sizeof(double));
    out->d = calloc(nchain*nout, sizeof(double));
    out->v = calloc(nchain*nout, sizeof(double));
    out->h = calloc(nchain*nout, sizeof(double));
    out->P = calloc(nchain*nout, sizeof(double));
    out->n = calloc(nchain*nout, sizeof(double));
    out->a = calloc(nchain*nout, sizeof(double));
    out->e = calloc(nchain*nout, sizeof(double));
    out->inc = calloc(nchain*nout, sizeof(double));
    out->Omega = calloc(nchain*nout, sizeof(double));
    out->omega = calloc(nchain*nout, sizeof(double));
    out->pomega = calloc(nchain*nout, sizeof(double));
    out->f = calloc(nchain*nout, sizeof(double));
    out->M = calloc(nchain*nout, sizeof(double));
    out->l = calloc(nchain*nout, sizeof(double));
    out->theta = calloc(nchain*nout, sizeof(double));
    out->T = calloc(nchain*nout, sizeof(double));
    out->rhill = calloc(nchain*nout, sizeof(double));
    out->Sigma = calloc(nchain*nout, sizeof(double));
}

void out_to_hdf5(struct output_structure* out){
    hid_t       file_id;
    //hsize_t     dims[RANK]={2,3};
    int datarank = 2;
    hsize_t scalardims[1] = {out->nout};
    hsize_t dims[2];
    dims[0] = out->nchain;
    dims[1] = out->nout;

    /* create a HDF5 file */
    file_id = H5Fcreate ("orbits.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t rootid = H5Gopen(file_id, "/", H5P_DEFAULT);
    H5LTset_attribute_int( rootid, "/", "nchain", &out->nchain, 1);
    H5LTset_attribute_int( rootid, "/", "nout",   &out->nout, 1);
    H5LTset_attribute_int( rootid, "/", "p",      &out->p, 1);
    H5LTset_attribute_int( rootid, "/", "q",      &out->q, 1);
    H5LTset_attribute_double( rootid, "/", "pmass", out->pmass, nchain);
    H5LTset_attribute_double( rootid, "/", "ai",  &out->ai, 1);


    H5LTset_attribute_double( rootid, "/", "Sigma0",   &out->Sigma0, 1);
    H5LTset_attribute_double( rootid, "/", "redge",   &out->redge, 1);
    H5LTset_attribute_double( rootid, "/", "deltaredge",   &out->deltaredge, 1);
    H5LTset_attribute_double( rootid, "/", "aspectratio0",   &out->aspectratio0, 1);
    H5LTset_attribute_double( rootid, "/", "alpha",   &out->alpha, 1);
    H5LTset_attribute_double( rootid, "/", "flaringindex",   &out->flaringindex, 1);
    H5LTset_attribute_double( rootid, "/", "ffudge",   &out->ffudge, 1);
    H5LTset_attribute_double( rootid, "/", "tmax",   &out->tmax, 1);
    H5LTset_attribute_double( rootid, "/", "tdep",   &out->tdep, 1);
    H5LTset_attribute_double( rootid, "/", "deltatdep",   &out->deltatdep, 1);

    H5LTmake_dataset_double(file_id, "/t", 1, scalardims, out->t);
    H5LTmake_dataset_double(file_id, "/d", datarank, dims, out->d);
    H5LTmake_dataset_double(file_id, "/v", datarank, dims, out->v);
    H5LTmake_dataset_double(file_id, "/h", datarank, dims, out->h);
    H5LTmake_dataset_double(file_id, "/P", datarank, dims, out->P);
    H5LTmake_dataset_double(file_id, "/n", datarank, dims, out->n);
    H5LTmake_dataset_double(file_id, "/a", datarank, dims, out->a);
    H5LTmake_dataset_double(file_id, "/e", datarank, dims, out->e);
    H5LTmake_dataset_double(file_id, "/inc", datarank, dims, out->inc);
    H5LTmake_dataset_double(file_id, "/Omega", datarank, dims, out->Omega);
    H5LTmake_dataset_double(file_id, "/omega", datarank, dims, out->omega);
    H5LTmake_dataset_double(file_id, "/pomega", datarank, dims, out->pomega);
    H5LTmake_dataset_double(file_id, "/f", datarank, dims, out->f);
    H5LTmake_dataset_double(file_id, "/M", datarank, dims, out->M);
    H5LTmake_dataset_double(file_id, "/l", datarank, dims, out->l);
    H5LTmake_dataset_double(file_id, "/theta", datarank, dims, out->theta);
    H5LTmake_dataset_double(file_id, "/T", datarank, dims, out->T);
    H5LTmake_dataset_double(file_id, "/rhill", datarank, dims, out->rhill);
    H5LTmake_dataset_double(file_id, "/Sigma", datarank, dims, out->Sigma);
   
    /* close file */
    H5Fclose (file_id);
}

void free_output_structure(struct output_structure* out){
   free(out->pmass);
   free(out->t);
   free(out->d);
   free(out->v);
   free(out->h);
   free(out->P);
   free(out->n);
   free(out->a);
   free(out->e);
   free(out->inc);
   free(out->Omega);
   free(out->omega);
   free(out->pomega);
   free(out->f);
   free(out->M);
   free(out->l);
   free(out->theta);
   free(out->T);
   free(out->rhill);
   free(out->Sigma);
}


double disc_surface_density(const double r, const double t){
    double  Sigma = out.Sigma0 * pow(r, -out.alpha) * tanh( (r - out.redge) / out.deltaredge) ; 
    if (t > tdep) {
        if (t < tdep + deltatdep){
            Sigma *= (1.0 - (t - out.tdep)/(out.deltatdep));
        } else {
            Sigma = 0.0;
        }
    } 
    return Sigma;
}

void migration_forces_tanaka_ward(struct reb_simulation* r){
    const double G = r->G;
    const int N = r->N;
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0]; // calculate migration forces with respect to center of mass;
    for(int i=1;i<N;i++){
        struct reb_particle* p = &(particles[i]);
        const double dvx = p->vx - com.vx;
        const double dvy = p->vy - com.vy;
        const double dvz = p->vz - com.vz;
        const double dx  = p->x  - com.x;
        const double dy  = p->y  - com.y;
        const double dz  = p->z  - com.z;
        const double rp = sqrt( dx*dx + dy*dy + dz*dz);
        const double theta = atan2(dy,dx);
        const double dvr     = dvx * (cos(theta))  + dvy * sin(theta);
        const double dvtheta = dvx * (-sin(theta)) + dvy * cos(theta);
        const double Omega = sqrt( G * com.m / pow(rp,3) );
        const double vk = Omega * rp;
        const double aspectratio = out.aspectratio0 * pow(rp, out.flaringindex);
        const double A_c_r     =  0.057; 
        const double A_c_theta = -0.868;
        const double A_c_z     = -1.088;
        const double A_s_r     =  0.176;
        const double A_s_theta =  0.325;
        const double A_s_z     = -0.871;

        const double Sigma = disc_surface_density(rp, r->t);

        // orignal formulae have vk/cs all over - replace with aspectratio = cs/vk
        const double fdampr     = (p->m/com.m) * pow(aspectratio,-4) *(Sigma *rp*rp / com.m) * Omega
                                  *(2.0 * A_c_r * (dvtheta - rp * Omega) + A_s_r * dvr);
        const double fdamptheta = (p->m/com.m) * pow(aspectratio,-4) *(Sigma *rp*rp / com.m) * Omega
                                  *(2.0 * A_c_theta * (dvtheta - rp * Omega) + A_s_theta * dvr);
        const double fdampz     = (p->m/com.m) * pow(aspectratio,-4) *(Sigma *rp*rp / com.m) * Omega
                                  *(2.0 * A_c_z * dvz + A_s_z * dz * Omega);
        const double fmigtheta  = out.ffudge * -2.17 * (p->m/com.m) * pow(aspectratio,-2) 
                                  * (Sigma *rp*rp / com.m) * Omega * vk;
        const double fdampx = cos(theta) * fdampr - sin(theta) * fdamptheta;
        const double fdampy = sin(theta) * fdampr + cos(theta) * fdamptheta;
        const double fmigx  = - sin(theta) * fmigtheta;
        const double fmigy  = + cos(theta) * fmigtheta;
        p-> ax += fdampx + fmigx;
        p-> ay += fdampy + fmigy;
        p-> az += fdampz;
    }
}

double get_particle_r(struct reb_particle* p, struct reb_particle* com){
    const double dx  = p->x -com->x;
    const double dy  = p->y -com->y;
    const double dz  = p->z -com->z;
    const double rp = sqrt( dx*dx + dy*dy + dz*dz);
    return rp;
}

void heartbeat(struct reb_simulation* r){

    if(reb_output_check(r, 20.*M_PI)){
        reb_output_timing(r, tmax);
    }
    if(reb_output_check(r, dtout)){
        reb_integrator_synchronize(r);
        reb_output_orbits(r,"orbits.txt");
        const double G = r->G;
        const int N = r->N;
        struct reb_particle* const particles = r->particles;
        struct reb_particle primary = particles[0];
        out.t[iout] = r->t;
        for(int i=1;i<N;i++){
            int errcode = 0;
            struct reb_orbit orb = reb_tools_particle_to_orbit_err(G, particles[i], primary, &errcode);

            out.d[(i-1)*out.nout +iout] = orb.d;
            out.v[(i-1)*out.nout +iout] = orb.v;
            out.h[(i-1)*out.nout +iout] = orb.h;
            out.P[(i-1)*out.nout +iout] = orb.P;
            out.n[(i-1)*out.nout +iout] = orb.n;
            out.a[(i-1)*out.nout +iout] = orb.a;
            out.e[(i-1)*out.nout +iout] = orb.e;
            out.inc[(i-1)*out.nout +iout] = orb.inc;
            out.Omega[(i-1)*out.nout +iout] = orb.Omega;
            out.omega[(i-1)*out.nout +iout] = orb.omega;
            out.pomega[(i-1)*out.nout +iout] = orb.pomega;
            out.f[(i-1)*out.nout +iout] = orb.f;
            out.M[(i-1)*out.nout +iout] = orb.M;
            out.l[(i-1)*out.nout +iout] = orb.l;
            out.theta[(i-1)*out.nout +iout] = orb.theta;
            out.T[(i-1)*out.nout +iout] = orb.T;
            out.rhill[(i-1)*out.nout +iout] = orb.rhill;
            const double rp = get_particle_r(&particles[i], &primary);
            out.Sigma[(i-1)*out.nout +iout] = disc_surface_density(rp, r->t);
        }
        iout += 1;
        reb_move_to_com(r); 
    }
    // should save simulation archives of snapshots here for a slected bunch of random times?
    //no, just randomize the evaporations scheduale
}
