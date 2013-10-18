#ifndef MSVMMAJ_TRAIN_DATASET_H
#define MSVMMAJ_TRAIN_DATASET_H

#include "globals.h"
#include "types.h"

struct Task {
	KernelType kerneltype;
	int weight_idx;
	long folds;
	long ID;
	double p;
	double kappa;
	double lambda;
	double epsilon;
	double *kernel_param;
	struct MajData *train_data;
	struct MajData *test_data;
	double performance;
};

struct Queue {
	struct Task **tasks;
	long N;
	long i;
};

struct Training {
	TrainType traintype;
	long repeats;
	long folds;
	long Np;
	long Nl;
	long Nk;
	long Ne;
	long Nw;
	int *weight_idxs;
	double *ps;
	double *lambdas;
	double *kappas;
	double *epsilons;
	char *train_data_file;
	char *test_data_file;
};

void make_queue(struct Training *training, struct Queue *queue,
		struct MajData *train_data, struct MajData *test_data);

struct Task *get_next_task(struct Queue *q);
void start_training_tt(struct Queue *q);
void start_training_cv(struct Queue *q);
void free_queue(struct Queue *q);

void consistency_repeats(struct Queue *q, long repeats, TrainType traintype);

double cross_validation(struct MajModel *model, struct MajModel *seed_model,
		struct MajData *data, long folds);

void make_model_from_task(struct Task *task, struct MajModel *model);
void copy_model(struct MajModel *from, struct MajModel *to);
#endif
