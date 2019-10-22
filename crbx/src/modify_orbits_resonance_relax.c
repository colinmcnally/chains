/**
 * @file    modify_orbits_forces.c
 * @brief   Update orbital elements with prescribed timescales using forces.
 * @author  Colin McNally <colin@colinmcnally.ca>
 * 
 * @section     LICENSE
 * Copyright (c) 2019 Colin McNally
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Orbit Modifications$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 C.P. McNally
 * Implementation Paper    *In progress*
 * Based on                `Matsumoto et al. 2012 2012Icar..221..624M`_.
 * C Example               
 * Python Example          
 * ======================= ===============================================
 * 
 * This applies physical forces that orbit-average to give exponential growth/decay of the semimajor axis, eccentricity and inclination.
 * 
 * **Effect Parameters**
 *
 * If coordinates not, defaults to using Jacobi coordinates.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * coordinates (enum)           No          Type of elements to use for modification (Jacobi, barycentric or particle).
 *                                          See the examples for usage.
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * One can pick and choose which particles have which parameters set.  
 * For each particle, any unset parameter is ignored.
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * res_aspectratio0 (double)    Yes         Disc aspect ratio at r=1
 * res_flaringindex (double)    Yes         Disc flaring index
 * res_ffudge (double)          Yes         Fudge factor controlling ratio of semimajor axis to eccentricity damping         
 * res_sigma0 (double)          Yes         Disc surface density at r=1         
 * res_alpha (double)           Yes         Disc negative surface density slope Sigma = Sigma0 r^-alpha   
 * res_redge (double)           Yes         Disc inner edge radius         
 * res_deltaredge (double)      Yes         Disc inner edge tanh transition width      
 * res_tdep (double)            Yes         Time at which the disc starts depleting
 * res_deltatdep (double)       Yes         Time over which disc depletion is finished
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

static struct reb_vec3d rebx_calculate_src_modify_orbits_resonance_relax(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    const double G = sim->G;
    double aspectratio0 = INFINITY;
    double flaringindex = INFINITY;
    double ffudge       = INFINITY;
    double sigma0       = INFINITY;
    double alpha        = INFINITY;
    double redge        = INFINITY;
    double deltaredge   = INFINITY;
    double tdep         = INFINITY;
    double deltatdep    = INFINITY;
 
    const double* const aspectratio0_ptr = rebx_get_param(sim->extras, force->ap, "res_aspectratio0");
    const double* const flaringindex_ptr = rebx_get_param(sim->extras, force->ap, "res_flaringindex");
    const double* const ffudge_ptr       = rebx_get_param(sim->extras, force->ap, "res_ffudge");
    const double* const sigma0_ptr       = rebx_get_param(sim->extras, force->ap, "res_sigma0");
    const double* const alpha_ptr        = rebx_get_param(sim->extras, force->ap, "res_alpha");
    const double* const redge_ptr        = rebx_get_param(sim->extras, force->ap, "res_redge");
    const double* const deltaredge_ptr   = rebx_get_param(sim->extras, force->ap, "res_deltaredge");
    const double* const tdep_ptr         = rebx_get_param(sim->extras, force->ap, "res_tdep");
    const double* const deltatdep_ptr    = rebx_get_param(sim->extras, force->ap, "res_deltatdep");

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x - source->x;
    const double dy = p->y - source->y;
    const double dz = p->z - source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;
    const double rp = sqrt(r2);

    struct reb_vec3d a = {0};

    if(aspectratio0_ptr != NULL){
        aspectratio0 = *aspectratio0_ptr;
    } else {
        reb_error(sim,"aspectratio0 not specified\n");
        return a;
    }

    if(flaringindex_ptr != NULL){
        flaringindex = *flaringindex_ptr;
    } else {
        reb_error(sim,"flaringindex not specified\n");
        return a;
    }

    if(ffudge_ptr != NULL){
        ffudge = *ffudge_ptr;
    } else {
        reb_error(sim,"ffudge not specified\n");
        return a;
    }

    if(sigma0_ptr != NULL){
        sigma0 = *sigma0_ptr;
    } else {
        reb_error(sim,"sigma0 not specified\n");
        return a;
    }

    if(alpha_ptr != NULL){
        alpha = *alpha_ptr;
    } else {
        reb_error(sim,"alpha not specified\n");
        return a;
    }

    if(redge_ptr != NULL){
        redge = *redge_ptr;
    } else {
        reb_error(sim,"redge not specified\n");
        return a;
    }

    if(deltaredge_ptr != NULL){
        deltaredge = *deltaredge_ptr;
    } else {
        reb_error(sim,"deltaredge not specified\n");
        return a;
    }

    if(tdep_ptr != NULL){
        tdep = *tdep_ptr;
    } else {
        reb_error(sim,"tdep not specified\n");
        return a;
    }

    if(deltatdep_ptr != NULL){
        deltatdep = *deltatdep_ptr;
    } else {
        reb_error(sim,"deltatdep not specified\n");
        return a;
    }


    double sigma = sigma0 * pow(rp, -alpha) * tanh( (rp - redge) / deltaredge); 
    if (sim->t > tdep){
        if (sim->t < tdep + deltatdep){
            sigma *= (1.0 - (sim->t - tdep)/deltatdep);
        } else {
            sigma = 0.0;
        }
    }

    const double theta = atan2(dy,dx);
    const double dvr     = dvx * (cos(theta))  + dvy * sin(theta);
    const double dvtheta = dvx * (-sin(theta)) + dvy * cos(theta);
    const double Omega = sqrt( G * source->m / pow(rp,3) );
    const double vk = Omega * rp;
    const double aspectratio = aspectratio0 * pow(rp, flaringindex);
    const double A_c_r     =  0.057; 
    const double A_c_theta = -0.868;
    const double A_c_z     = -1.088;
    const double A_s_r     =  0.176;
    const double A_s_theta =  0.325;
    const double A_s_z     = -0.871;

    // orignal formulae have vk/cs all over - replace with aspectratio = cs/vk
    const double fdampr     = (p->m/source->m) * pow(aspectratio,-4) *(sigma *rp*rp / source->m) * Omega
                              *(2.0 * A_c_r * (dvtheta - rp * Omega) + A_s_r * dvr);
    const double fdamptheta = (p->m/source->m) * pow(aspectratio,-4) *(sigma *rp*rp / source->m) * Omega
                              *(2.0 * A_c_theta * (dvtheta - rp * Omega) + A_s_theta * dvr);
    const double fdampz     = (p->m/source->m) * pow(aspectratio,-4) *(sigma *rp*rp / source->m) * Omega
                              *(2.0 * A_c_z * dvz + A_s_z * dz * Omega);
    const double fmigtheta  = ffudge * -2.17 * (p->m/source->m) * pow(aspectratio,-2) 
                              * (sigma *rp*rp / source->m) * Omega * vk;
    const double fdampx = cos(theta) * fdampr - sin(theta) * fdamptheta;
    const double fdampy = sin(theta) * fdampr + cos(theta) * fdamptheta;
    const double fmigx  = - sin(theta) * fmigtheta;
    const double fmigy  = + cos(theta) * fmigtheta;

    a.x += fdampx + fmigx;
    a.y += fdampy + fmigy;
    a.z += fdampz;

    return a;
}

void rebx_modify_orbits_resonance_relax(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_src_modify_orbits_resonance_relax, particles, N);
}
