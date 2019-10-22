/**
 * @file    modify_orbits_forces_edge.c
 * @brief   Update orbital elements with prescribed timescales using forces, inclide innder edge.
 * @author  Colin McNally <colin@colinmcnally.ca>, Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
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
 * Authors                 C. McNally from original of D. Tamayo, H. Rein
 * Implementation Paper    *In progress*
 * Based on                `Papaloizou & Larwood 2000 <http://labs.adsabs.harvard.edu/adsabs/abs/2000MNRAS.315..823P/>`_.
 * C Example               
 * Python Example          
 * ======================= ===============================================
 * 
 * This applies physical forces that orbit-average to give exponential growth/decay of the semimajor axis, eccentricity and inclination.
 * The eccentricity damping keeps the angular momentum constant (corresponding to `p=1` in modify_orbits_direct), which means that eccentricity damping will induce some semimajor axis evolution.
 * Additionally, eccentricity/inclination damping will induce pericenter/nodal precession.
 * Both these effects are physical, and the method is more robust for strongly perturbed systems.
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
 * tau_a (double)               No          Semimajor axis exponential growth/damping timescale
 * tau_e (double)               No          Eccentricity exponential growth/damping timescale
 * tau_inc (double)             No          Inclination axis exponential growth/damping timescale
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

static struct reb_vec3d rebx_calculate_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source){
    double tau_a        = INFINITY;
    double tau_e        = INFINITY;
    double tau_inc      = INFINITY;
    double redge        = INFINITY;
    double deltaredge   = INFINITY;
    double tdep         = INFINITY;
    double deltatdep    = INFINITY;
    
    const double* const tau_a_ptr        = rebx_get_param(sim->extras, p->ap, "tau_a");
    const double* const tau_e_ptr        = rebx_get_param(sim->extras, p->ap, "tau_e");
    const double* const tau_inc_ptr      = rebx_get_param(sim->extras, p->ap, "tau_inc");
    const double* const redge_ptr        = rebx_get_param(sim->extras, force->ap, "res_redge");
    const double* const deltaredge_ptr   = rebx_get_param(sim->extras, force->ap, "res_deltaredge");
    const double* const tdep_ptr         = rebx_get_param(sim->extras, force->ap, "res_tdep");
    const double* const deltatdep_ptr    = rebx_get_param(sim->extras, force->ap, "res_deltatdep");

    const double dvx = p->vx - source->vx;
    const double dvy = p->vy - source->vy;
    const double dvz = p->vz - source->vz;
    const double dx = p->x-source->x;
    const double dy = p->y-source->y;
    const double dz = p->z-source->z;
    const double r2 = dx*dx + dy*dy + dz*dz;
    const double rp = sqrt(r2);
    
    struct reb_vec3d a = {0};

    if(tau_a_ptr != NULL){
        tau_a = *tau_a_ptr;
    }
    if(tau_e_ptr != NULL){
        tau_e = *tau_e_ptr;
    }
    if(tau_inc_ptr != NULL){
        tau_inc = *tau_inc_ptr;
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

    double sigma = 1.0 - exp(-(rp-redge)/deltaredge); 
    if (sim->t > tdep){
        if (sim->t < tdep + deltatdep){
            sigma *= (1.0 - (sim->t - tdep)/deltatdep);
        } else {
            sigma = 0.0;
        }
    }

    a.x =  sigma*dvx/(2.*tau_a);
    a.y =  sigma*dvy/(2.*tau_a);
    a.z =  sigma*dvz/(2.*tau_a);

    if (tau_e < INFINITY || tau_inc < INFINITY){
        const double vdotr = dx*dvx + dy*dvy + dz*dvz;
        const double prefac = 2*vdotr/r2/tau_e;
        a.x += prefac*dx;
        a.y += prefac*dy;
        a.z += prefac*dz + 2.*dvz/tau_inc;
    }
    return a;
}

void rebx_modify_orbits_forces_edge(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    int* ptr = rebx_get_param(sim->extras, force->ap, "coordinates");
    enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI; // Default
    if (ptr != NULL){
        coordinates = *ptr;
    }
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebx_com_force(sim, force, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_forces, particles, N);
}
