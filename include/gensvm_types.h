/**
 * @file gensvm_types.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Definitions of common types
 *
 * @details
 * Here common types used throughout the program are defined.
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

#ifndef GENSVM_TYPES_H
#define GENSVM_TYPES_H

/**
 * @brief type of training used in parameter grid search
 */
typedef enum {
	CV=0, /**< cross validation */
	TT=1  /**< data with existing train/test split */
} TrainType;

/**
 * @brief type of kernel used in training
 */
typedef enum {
	K_LINEAR=0, 	/**< Linear kernel */
	K_POLY=1, 	/**< Polynomial kernel */
	K_RBF=2, 	/**< RBF kernel */
	K_SIGMOID=3,  	/**< Sigmoid kernel */
} KernelType;

#endif
