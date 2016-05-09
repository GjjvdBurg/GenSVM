/**
 * @file types.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Definitions of common types
 *
 * @details
 * Here common types used throughout the program are defined.
 *
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
