/**
 * @file gensvm_timer.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Header file for gensvm_timer.c
 *
 * @details
 * Function declaration for timer function used to measure computation time.
 *
 * @copyright
 Copyright 2016, G.J.J. van den Burg.

 This file is part of GenSVM.

 GenSVM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GenSVM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GenSVM. If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef GENSVM_TIMER_H
#define GENSVM_TIMER_H

// includes
#include "gensvm_globals.h"

// macro for the GenTimer type
#ifdef ON_WINDOWS
#define GenTime LARGE_INTEGER
#else
#define GenTime struct timeval
#endif

/// Timer macro for easily recording time
#define Timer(spec) gensvm_set_time(&spec)

// function declarations
void gensvm_set_time(GenTime *t);
double gensvm_elapsed_time(GenTime *start, GenTime *stop);
void gensvm_sleep(double seconds);

#endif
