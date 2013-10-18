#ifndef TYPES_H
#define TYPES_H

typedef enum {
	false,
	true
} bool;

typedef enum {
	CV=0,
	TT=1
} TrainType;

typedef enum {
	K_LINEAR=0,
	K_POLY=1,
	K_RBF=2,
	K_SIGMOID=3,
} KernelType;

#endif
