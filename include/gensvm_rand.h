/**
 * @file gensvm_rand.h
 * @author G.J.J. van den Burg
 * @date 2018-04-17
 * @brief Header files for gensvm_rand.c
 *
 * @copyright
 * Copyright 2018, G.J.J. van den Burg

 GenSVM is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 GenSVM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GenSVM; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#ifndef GENSVM_RAND_H
#define GENSVM_RAND_H

#include "gensvm_globals.h"

// structs
struct gensvm_random_data
{
	int32_t *fptr;		/* Front pointer.  */
	int32_t *rptr;		/* Rear pointer.  */
	int32_t *state;		/* Array of state values.  */
	int rand_type;		/* Type of random number generator.  */
	int rand_deg;		/* Degree of random number generator.  */
	int rand_sep;		/* Distance between front and rear.  */
	int32_t *end_ptr;	/* Pointer behind state table.  */
};

// forward declarations
int gensvm_rand(void);
void gensvm_srand(unsigned int x);

#endif
