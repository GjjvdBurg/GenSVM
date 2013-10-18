#ifndef MSVMMAJ_TRAIN_H
#define MSVMMAJ_TRAIN_H

#include "globals.h"

//forward declarations
struct MajData;
struct MajModel;

// function declarations
void msvmmaj_optimize(struct MajModel *model, struct MajData *data);

double msvmmaj_get_loss(struct MajModel *model, struct MajData *data,
		double *ZV);

void msvmmaj_get_update(struct MajModel *model, struct MajData *data,
		double *B, double *ZAZ, double *ZAZV, double *ZAZVT);

#endif
