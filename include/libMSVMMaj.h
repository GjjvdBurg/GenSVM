#ifndef LIBMSVMMAJ_H
#define LIBMSVMMAJ_H

#include "globals.h"

// forward declarations 
struct MajData;
struct MajModel;

// function declarations
void msvmmaj_simplex_gen(long K, double *U);
void msvmmaj_category_matrix(struct MajModel *model, struct MajData *data);
void msvmmaj_simplex_diff(struct MajModel *model, struct MajData *dataset);

void msvmmaj_calculate_errors(struct MajModel *model, struct MajData *data, double *ZV);
void msvmmaj_calculate_huber(struct MajModel *model);

void msvmmaj_step_doubling(struct MajModel *model);

void msvmmaj_seed_model_V(struct MajModel *from_model, struct MajModel *to_model);
void msvmmaj_initialize_weights(struct MajData *data, struct MajModel *model);

#endif
