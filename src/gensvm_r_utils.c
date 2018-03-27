/**
 * @file gensvm_r_utils.c
 * @author G.J.J. van den Burg
 * @date 2018-03-26
 * @brief Utilities for dealing with interrupts from R cleanly

 * Copyright (C) G.J.J. van den Burg

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

 */

#define STRICT_R_HEADERS

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
// from: https://stackoverflow.com/q/40563522

#include "gensvm_r_utils.h"

bool __R_requested_interrupt = false;

void _interrupt_fn(void *dummy)
{
	R_CheckUserInterrupt();
}

void gensvm_R_reset_interrupt_hdl(void)
{
	__R_requested_interrupt = false;
}

bool gensvm_R_pending_interrupt(void)
{
	if (__R_requested_interrupt)
		return true;
	__R_requested_interrupt = R_ToplevelExec(_interrupt_fn, NULL) == FALSE;
	return __R_requested_interrupt;
}
