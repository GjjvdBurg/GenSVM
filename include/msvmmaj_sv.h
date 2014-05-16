/**
 * @file msvmmaj_sv.h
 * @author Gertjan van den Burg
 * @date May, 2014
 * @brief Header file for msvmmaj_sv.c
 *
 * @details
 * Contains function declarations for functions used to count support vectors.
 *
 */

#ifndef MSVMMAJ_SV_H
#define MSVMMAJ_SV_H

#include "globals.h"

long msvmmaj_num_sv(struct MajModel *model, struct MajData *data);

#endif
