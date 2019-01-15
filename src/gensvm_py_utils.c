/**
 * @file gensvm_py_utils.c
 * @author G.J.J. van den Burg
 * @date 2018-04-30
 * @brief Utilities for dealing with interrupts from Python cleanly

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

// from: https://stackoverflow.com/q/40563522

#include "gensvm_py_utils.h"

bool __py_requested_interrupt = false;

void gensvm_py_reset_interrupt_hdl(void)
{
	__py_requested_interrupt = false;
}

bool gensvm_py_pending_interrupt(void)
{
	if (__py_requested_interrupt)
		return true;

	PyObject* pyerr = PyErr_Occurred();
	__py_requested_interrupt = (pyerr != NULL);
	return __py_requested_interrupt;
}
