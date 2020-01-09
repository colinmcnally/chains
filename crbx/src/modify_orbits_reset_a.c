/**
 * @file    modify_orbits_reset_a.c
 * @brief   Update orbital set semimajor axis at each timestep..
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2020 Colin McNally
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
 * Authors                 C.McNally, from original by D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                `Lee & Peale 2002 <http://labs.adsabs.harvard.edu/adsabs/abs/2002ApJ...567..596L/>`_. 
 * C Example               
 * Python Example          
 * ======================= ===============================================
 * 
 * Aritifically reset the semimajor axis of a particle.
 * 
 * **Effect Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * set_a (double)               No          Semimajor axis to set particle to
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static struct reb_particle rebx_calculate_modify_orbits_reset_a(struct reb_simulation* const sim, struct rebx_operator* const operator, struct reb_particle* p, struct reb_particle* primary, const double dt){
    struct rebx_extras* const rebx = sim->extras;
    int err=0;
    double tdep         = INFINITY;
    double deltatdep    = INFINITY;
    struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *primary, &err);
    if(err){        // mass of primary was 0 or p = primary.  Return same particle without doing anything.
        return *p;
    }
    const double* const set_a            = rebx_get_param(rebx, p->ap, "set_a");
    const double* const tdep_ptr         = rebx_get_param(rebx, p->ap, "res_tdep");
    const double* const deltatdep_ptr    = rebx_get_param(rebx, p->ap, "res_deltatdep");
    
    if(set_a != NULL){

        if(tdep_ptr != NULL){
            tdep = *tdep_ptr;
        } else {
            reb_error(sim,"tdep not specified for a reset_a particle\n");
            return *p;
        }
        if(deltatdep_ptr != NULL){
            deltatdep = *deltatdep_ptr;
        } else {
            reb_error(sim,"deltatdep not specified for a reset_a particle\n");
            return *p;
        }
        if (sim->t < tdep + deltatdep)
            o.a = *set_a;
    }
    return reb_tools_orbit_to_particle(sim->G, *primary, p->m, o.a, o.e, o.inc, o.Omega, o.omega, o.f);
}

void rebx_modify_orbits_reset_a(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    const int* const ptr = rebx_get_param(sim->extras, operator->ap, "coordinates");
   	enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI;
	if (ptr != NULL){
		coordinates = *ptr;
	}
    const int back_reactions_inclusive = 1;
    const char* reference_name = "primary";
    rebxtools_com_ptm(sim, operator, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_reset_a, dt);
}
