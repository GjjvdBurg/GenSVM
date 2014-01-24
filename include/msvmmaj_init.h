/**
 * @file msvmmaj_init.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header file for msvmmaj_init.c
 *
 * @details
 * Contains function declarations for the initialization functions for
 * MajModel and MajData structures.
 */

#ifndef MSVMMAJ_INIT_H
#define MSVMMAJ_INIT_H

// forward declaration
struct MajData;
struct MajModel;

struct MajModel *msvmmaj_init_model();

struct MajData *msvmmaj_init_data();

void msvmmaj_allocate_model(struct MajModel *model);
void msvmmaj_free_model(struct MajModel *model);
void msvmmaj_free_data(struct MajData *data);

#endif
