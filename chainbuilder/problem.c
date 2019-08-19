// For generating and evolving resonant chain systems
// Derived from the RREBOUND example for migrating planets
// Keeps pretty close to the methods of Yuji Matsumoto, Makiko Nagasawa, Shigeru Ida 
//                                      2012Icar..221..624M
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>
#include "rebound.h"
#include "hdf5.h"
#include "hdf5_hl.h"

// snapshot interval in seconds
#define SNAP_WALLTIME_INTERVAL (60*30)

struct output_structure {
  int nout;
  int iout;
  int nchain;
  int p;
  int q;
  char seqstr[256];
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
  double dtout;
  double tdep;
  double deltatdep;
  double tcollstop; ///<Init to -1.0 and set if stopped due to collision
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

double next_snap_walltime;

bool file_exists(const char* file);
int collision_stop(struct reb_simulation* const r, struct reb_collision c);
void migration_forces(struct reb_simulation* r);
void migration_forces_tanaka_ward(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);
void init_output_structure(struct output_structure* out, const int nchain, const int p, const double tmax,
             const double tdep, const double deltatdep, const int seqnum);
void out_to_hdf5(struct output_structure* out);
void out_from_hdf5(struct output_structure* out);
void output_snapshot(struct output_structure* const out, struct reb_simulation* const r);
void free_output_structure(struct output_structure* out);
double disc_surface_density(const double r, const double t);
double get_particle_r(struct reb_particle* p, struct reb_particle* com);
void run_sim(const int nchain, const int p, const double tmax,
             const double tdep, const double deltatdep, const int seqnum);

int main(int argc, char* argv[]){
    // set output counter
    if (argc != 7){
        printf("Wrong number of arguments %i should be 7\n", argc);
        exit(-1);
    }
    run_sim( atoi(argv[1]), atoi(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atoi(argv[6]));
}

void run_sim(const int nchain, const int p, const double tmax,
             const double tdep, const double deltatdep, const int seqnum){
    init_output_structure(&out, nchain, p, tmax, tdep, deltatdep, seqnum);

    struct reb_simulation* r = reb_create_simulation();
    // Initial conditions
    // Parameters are those of Lee & Peale 2002, Figure 4. 
    struct reb_particle star = {0};
    star.m  = out.starmass;            // This is a sub-solar mass star
    reb_add(r, star); 
    
    double ai = out.ai;
    int errcode;
    for (int ip=0; ip < out.nchain; ip++){
        // struct reb_particle reb_tools_orbit_to_particle_err(double G, struct reb_particle primary, 
        //      double m, double a, double e, double i, double Omega, double omega, double f, int* err)
        //randomize inclination, eccentricity, true anomaly - from seed
        struct reb_particle pt = reb_tools_orbit_to_particle_err( 
                r->G, star, out.pmass[ip], ai, 0.01, 0.01, 0.0, 0.0, (double)ip, &errcode);
        double rhill = ai * pow(pt.m/(3.*star.m),1./3.);    // Hill radius
        pt.r = rhill; // Set planet radius to hill radius
                      // A collision is recorded when planets get within their hill radius
                      // The hill radius of the particles might change, so it should be recalculated after a while
        pt.lastcollision = 0;
        reb_add(r, pt);
        ai = ai * pow(out.p / (out.p + out.q + 0.05), -2./3.);
    } 


    // if a snapshot exists, just do a restart
    char filename[512];
    sprintf(filename,"orbits_%i_%i:%i_%.2e_%s.snap", out.nchain, out.p+out.q, out.p, out.pmass[0], out.seqstr);
    if (file_exists(filename)) {
        printf("Restart triggered by %s\n", filename);
        r = reb_create_simulation_from_binary(filename);
        out_from_hdf5(&out);
        // now reset out.iout
        out.iout = (int)(r->t/out.dtout)+1;
        printf(" r->t %e out.tcollstop %e\n", r->t, out.tcollstop);
    }

    // Do all the function pointer stuff _after_ reload from snapshot!
    r->integrator    = REB_INTEGRATOR_WHFAST;
    //r->integrator    = REB_INTEGRATOR_IAS15;
    r->dt            = 1e-2*2.*M_PI *pow(0.1, 1.5);        // in year/(2*pi)
    r->collision        = REB_COLLISION_DIRECT;
    r->collision_resolve     = collision_stop;        // Set function pointer for collision recording.
    r->additional_forces = migration_forces_tanaka_ward;     //Set function pointer to add dissipative forces.
    r->heartbeat = heartbeat;  
    r->force_is_velocity_dependent = 1;

    next_snap_walltime = SNAP_WALLTIME_INTERVAL;

    // only actually run if not at tmax, or if tcollstop is not yet set to > 0
    if (out.tmax > r->t && out.tcollstop < 0.0 ){
        reb_move_to_com(r);          
        reb_integrate(r, out.tmax);
        output_snapshot(&out, r);
        out_to_hdf5(&out);
    }else{
        printf("This snapshot was already past the tmax of the simulation, or a collision has happened so exiting.\n");
    }
    free_output_structure(&out);

    printf("\n");
}

/**
 * Check if a file exists
 * @return true if and only if the file exists, false else
 */
bool file_exists(const char* file) {
    struct stat buf;
    return (stat(file, &buf) == 0);
}

void init_output_structure(struct output_structure* out,
                           const int nchain, const int p, const double tmax,
                           const double tdep, const double deltatdep, const int seqnum){

    // disc model
    const double Sigma0       = 3.8e-4 *0.5*1.4;
    const double redge        = 0.1;
    const double deltaredge   = 0.001;
    const double aspectratio0 = 0.035;
    const double alpha        = 1.0;
    const double flaringindex = 2.0/7.0;
    const double ffudge       = 1.0/100.0; //fudge factor for Matsumoto et al. 2012 - not sure on sign
    const double dtout        = 500.0 * 2.0*M_PI*pow(0.1,1.5);

    const double starmass = 1.0;

    const int q = 1;
    const double pmass = 1e-5;

    sprintf(out->seqstr, "%i", seqnum);
    
    double ai = 0.11;
    int nout = (int)((tmax+dtout)/dtout) +2;
    out->nout = nout;
    out->iout = 0;
    out->nchain = nchain;
    out->tmax = tmax;
    out->dtout = dtout;
    out->tdep = tdep;
    out->tcollstop = -1.0;
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
    hid_t       fileid;
    //hsize_t     dims[RANK]={2,3};
    int datarank = 2;
    hsize_t scalardims[1] = {out->nout};
    hsize_t dims[2];
    dims[0] = out->nchain;
    dims[1] = out->nout;

    char filename[500];
    sprintf(filename,"orbits_%i_%i:%i_%.2e_%s.h5", out->nchain, out->p+out->q, out->p, out->pmass[0], out->seqstr);
    /* create a HDF5 file */
    fileid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t rootid = H5Gopen(fileid, "/", H5P_DEFAULT);
    H5LTset_attribute_int( rootid, "/", "nchain", &out->nchain, 1);
    H5LTset_attribute_int( rootid, "/", "nout",   &out->nout, 1);
    H5LTset_attribute_int( rootid, "/", "lastout",&out->iout, 1);
    H5LTset_attribute_int( rootid, "/", "p",      &out->p, 1);
    H5LTset_attribute_int( rootid, "/", "q",      &out->q, 1);
    H5LTset_attribute_double( rootid, "/", "pmass", out->pmass, out->nchain);
    H5LTset_attribute_double( rootid, "/", "ai",  &out->ai, 1);


    H5LTset_attribute_double( rootid, "/", "Sigma0",   &out->Sigma0, 1);
    H5LTset_attribute_double( rootid, "/", "redge",   &out->redge, 1);
    H5LTset_attribute_double( rootid, "/", "deltaredge",   &out->deltaredge, 1);
    H5LTset_attribute_double( rootid, "/", "aspectratio0",   &out->aspectratio0, 1);
    H5LTset_attribute_double( rootid, "/", "alpha",   &out->alpha, 1);
    H5LTset_attribute_double( rootid, "/", "flaringindex",   &out->flaringindex, 1);
    H5LTset_attribute_double( rootid, "/", "ffudge",   &out->ffudge, 1);
    H5LTset_attribute_double( rootid, "/", "tmax",   &out->tmax, 1);
    H5LTset_attribute_double( rootid, "/", "dtout",   &out->dtout, 1);
    H5LTset_attribute_double( rootid, "/", "tdep",   &out->tdep, 1);
    H5LTset_attribute_double( rootid, "/", "deltatdep",   &out->deltatdep, 1);
    H5LTset_attribute_double( rootid, "/", "tcollstop",   &out->tcollstop, 1);

    H5Gclose (rootid);

    H5LTmake_dataset_double(fileid, "/t", 1, scalardims, out->t);
    H5LTmake_dataset_double(fileid, "/d", datarank, dims, out->d);
    H5LTmake_dataset_double(fileid, "/v", datarank, dims, out->v);
    H5LTmake_dataset_double(fileid, "/h", datarank, dims, out->h);
    H5LTmake_dataset_double(fileid, "/P", datarank, dims, out->P);
    H5LTmake_dataset_double(fileid, "/n", datarank, dims, out->n);
    H5LTmake_dataset_double(fileid, "/a", datarank, dims, out->a);
    H5LTmake_dataset_double(fileid, "/e", datarank, dims, out->e);
    H5LTmake_dataset_double(fileid, "/inc", datarank, dims, out->inc);
    H5LTmake_dataset_double(fileid, "/Omega", datarank, dims, out->Omega);
    H5LTmake_dataset_double(fileid, "/omega", datarank, dims, out->omega);
    H5LTmake_dataset_double(fileid, "/pomega", datarank, dims, out->pomega);
    H5LTmake_dataset_double(fileid, "/f", datarank, dims, out->f);
    H5LTmake_dataset_double(fileid, "/M", datarank, dims, out->M);
    H5LTmake_dataset_double(fileid, "/l", datarank, dims, out->l);
    H5LTmake_dataset_double(fileid, "/theta", datarank, dims, out->theta);
    H5LTmake_dataset_double(fileid, "/T", datarank, dims, out->T);
    H5LTmake_dataset_double(fileid, "/rhill", datarank, dims, out->rhill);
    H5LTmake_dataset_double(fileid, "/Sigma", datarank, dims, out->Sigma);
   
    /* close file */
    H5Fclose (fileid);
}

//This function is only for re-populating on restart, so assume the out structure is already initialized
void out_from_hdf5(struct output_structure* out){
    hid_t       fileid;
    //hsize_t     dims[RANK]={2,3};
    char filename[512];
    sprintf(filename,"orbits_%i_%i:%i_%.2e_%s.h5", out->nchain, out->p+out->q, out->p, out->pmass[0], out->seqstr);

    /* create a HDF5 file */
    fileid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t rootid = H5Gopen(fileid, "/", H5P_DEFAULT);
    H5LTget_attribute_int( rootid, "/", "nchain", &out->nchain);
    H5LTget_attribute_int( rootid, "/", "nout",   &out->nout);
    H5LTget_attribute_int( rootid, "/", "lastout",&out->iout);
    H5LTget_attribute_int( rootid, "/", "p",      &out->p);
    H5LTget_attribute_int( rootid, "/", "q",      &out->q);
    H5LTget_attribute_double( rootid, "/", "pmass", out->pmass);
    H5LTget_attribute_double( rootid, "/", "ai",  &out->ai);

    H5LTget_attribute_double( rootid, "/", "Sigma0",   &out->Sigma0);
    H5LTget_attribute_double( rootid, "/", "redge",   &out->redge);
    H5LTget_attribute_double( rootid, "/", "deltaredge",   &out->deltaredge);
    H5LTget_attribute_double( rootid, "/", "aspectratio0",   &out->aspectratio0);
    H5LTget_attribute_double( rootid, "/", "alpha",   &out->alpha);
    H5LTget_attribute_double( rootid, "/", "flaringindex",   &out->flaringindex);
    H5LTget_attribute_double( rootid, "/", "ffudge",   &out->ffudge);
    H5LTget_attribute_double( rootid, "/", "tmax",   &out->tmax);
    H5LTget_attribute_double( rootid, "/", "dtout",   &out->dtout);
    H5LTget_attribute_double( rootid, "/", "tdep",   &out->tdep);
    H5LTget_attribute_double( rootid, "/", "deltatdep",   &out->deltatdep);
    H5LTget_attribute_double( rootid, "/", "tcollstop",   &out->tcollstop);

    H5Gclose (rootid);

    H5LTread_dataset_double(fileid, "/t", out->t);
    H5LTread_dataset_double(fileid, "/d", out->d);
    H5LTread_dataset_double(fileid, "/v", out->v);
    H5LTread_dataset_double(fileid, "/h", out->h);
    H5LTread_dataset_double(fileid, "/P", out->P);
    H5LTread_dataset_double(fileid, "/n", out->n);
    H5LTread_dataset_double(fileid, "/a", out->a);
    H5LTread_dataset_double(fileid, "/e", out->e);
    H5LTread_dataset_double(fileid, "/inc", out->inc);
    H5LTread_dataset_double(fileid, "/Omega", out->Omega);
    H5LTread_dataset_double(fileid, "/omega", out->omega);
    H5LTread_dataset_double(fileid, "/pomega", out->pomega);
    H5LTread_dataset_double(fileid, "/f", out->f);
    H5LTread_dataset_double(fileid, "/M", out->M);
    H5LTread_dataset_double(fileid, "/l", out->l);
    H5LTread_dataset_double(fileid, "/theta", out->theta);
    H5LTread_dataset_double(fileid, "/T", out->T);
    H5LTread_dataset_double(fileid, "/rhill", out->rhill);
    H5LTread_dataset_double(fileid, "/Sigma", out->Sigma);
   
    /* close file */
    H5Fclose (fileid);
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

void output_snapshot(struct output_structure* const out, struct reb_simulation* const r){
    char filename[512];
    sprintf(filename,"orbits_%i_%i:%i_%.2e_%s.snap", out->nchain, out->p+out->q, out->p, out->pmass[0], out->seqstr);
    reb_output_binary(r, filename);
}

// Define our own collision resolve function.
int collision_stop(struct reb_simulation* const r, struct reb_collision c){
    struct reb_particle com = r->particles[0]; 
    struct reb_particle p1  = r->particles[c.p1];
    struct reb_particle p2  = r->particles[c.p2];
    
    double dx  = p1.x  - com.x;
    double dy  = p1.y  - com.y;
    double dz  = p1.z  - com.z;
    const double rp1 = sqrt( dx*dx + dy*dy + dz*dz);    
    const double rhill1 =  rp1 * pow(p1.m/(3.*com.m),1./3.);

    dx  = p2.x  - com.x;
    dy  = p2.y  - com.y;
    dz  = p2.z  - com.z;
    const double rp2 = sqrt( dx*dx + dy*dy + dz*dz);    
    const double rhill2 =  rp2 * pow(p2.m/(3.*com.m),1./3.);

    const double cx = p1.x - p2.x;
    const double cy = p1.y - p2.y;
    const double cz = p1.z - p2.z;
    const double rcoll = sqrt(cx*cx + cy*cy + cz*cz);

    // check if we are really within the current rhill at present position
    // this could be several other realted ccriteria, like mutual Hill radius
    if (rhill1 > rcoll || rhill2 > rcoll){
        r->status = REB_EXIT_COLLISION;
        r->particles[c.p1].lastcollision = r->t;
        r->particles[c.p2].lastcollision = r->t;
        out.tcollstop = r->t;
        printf("Stop on collision %e\n",r->t);
    }
    return 0; // don't remove either particle
}


double disc_surface_density(const double r, const double t){
    double  Sigma = out.Sigma0 * pow(r, -out.alpha) * tanh( (r - out.redge) / out.deltaredge) ; 
    if (t > out.tdep) {
        if (t < out.tdep + out.deltatdep){
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

    //r->walltime is a double of seconds since start in the step routine only
    if(r->walltime >= next_snap_walltime){
        output_snapshot(&out, r);
        out_to_hdf5(&out);
        next_snap_walltime += SNAP_WALLTIME_INTERVAL;
    }
    if(reb_output_check(r, 200.*M_PI)){
        reb_output_timing(r, out.tmax);
        printf("\n");
        printf(" walltime %e next snap %e\n",r->walltime, next_snap_walltime);
    }
    if(reb_output_check(r, out.dtout) && (r->t > 0.0 && r->t >= out.dtout)){
        reb_integrator_synchronize(r);
        //reb_output_orbits(r,"orbits.txt");
        const double G = r->G;
        const int N = r->N;
        struct reb_particle* const particles = r->particles;
        struct reb_particle primary = particles[0];
        out.t[out.iout] = r->t;
        assert(out.iout < out.nout);
        for(int i=1;i<N;i++){
            int errcode = 0;
            struct reb_orbit orb = reb_tools_particle_to_orbit_err(G, particles[i], primary, &errcode);

            out.d[(i-1)*out.nout +out.iout] = orb.d;
            out.v[(i-1)*out.nout +out.iout] = orb.v;
            out.h[(i-1)*out.nout +out.iout] = orb.h;
            out.P[(i-1)*out.nout +out.iout] = orb.P;
            out.n[(i-1)*out.nout +out.iout] = orb.n;
            out.a[(i-1)*out.nout +out.iout] = orb.a;
            out.e[(i-1)*out.nout +out.iout] = orb.e;
            out.inc[(i-1)*out.nout +out.iout] = orb.inc;
            out.Omega[(i-1)*out.nout +out.iout] = orb.Omega;
            out.omega[(i-1)*out.nout +out.iout] = orb.omega;
            out.pomega[(i-1)*out.nout +out.iout] = orb.pomega;
            out.f[(i-1)*out.nout +out.iout] = orb.f;
            out.M[(i-1)*out.nout +out.iout] = orb.M;
            out.l[(i-1)*out.nout +out.iout] = orb.l;
            out.theta[(i-1)*out.nout +out.iout] = orb.theta;
            out.T[(i-1)*out.nout +out.iout] = orb.T;
            out.rhill[(i-1)*out.nout +out.iout] = orb.rhill;
            const double rp = get_particle_r(&particles[i], &primary);
            out.Sigma[(i-1)*out.nout +out.iout] = disc_surface_density(rp, r->t);
        }
        //printf("%e iout %i nout %i %f %f %e %e\n",r->t, iout, out.nout, floor((r->t+out.dtout)/out.dtout), floor((r->t+out.dtout -r->dt)/out.dtout), r->t, r->dt);
        out.iout += 1;
        reb_move_to_com(r); 
    }
}
