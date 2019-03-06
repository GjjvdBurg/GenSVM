/**
 * @file gensvm_r_utils.h
 * @author G.J.J. van den Burg
 * @date 2018-03-26
 * @brief Header definitions for gensvm_r_utils.c

 * Copyright (C) 2018, G.J.J. van den Burg

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


#ifndef GENSVM_PY_UTILS_H
#define GENSVM_PY_UTILS_H

#include "gensvm_globals.h"

void gensvm_py_reset_interrupt_hdl(void);
bool gensvm_py_pending_interrupt(void);

#endif
