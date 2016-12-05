/**
 * @file test_gensvm_gridsearch.c
 * @author G.J.J. van den Burg
 * @date 2016-10-18
 * @brief Unit tests for gensvm_gridsearch.c
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

#include "minunit.h"
#include "gensvm_gridsearch.h"

extern FILE *GENSVM_OUTPUT_FILE;

char *test_fill_queue_nokernel()
{
	struct GenData *train_data = gensvm_init_data();
	struct GenData *test_data = gensvm_init_data();
	struct GenGrid *grid = gensvm_init_grid();

	grid->Np = 3;
	grid->Nl = 2;
	grid->Nk = 1;
	grid->Ne = 1;
	grid->Nw = 1;

	grid->ps = Calloc(double, grid->Np);
	grid->ps[0] = 1.0;
	grid->ps[1] = 1.5;
	grid->ps[2] = 2.0;

	grid->lambdas = Calloc(double, grid->Nl);
	grid->lambdas[0] = 1.0;
	grid->lambdas[1] = 5.0;

	grid->kappas = Calloc(double, grid->Nk);
	grid->kappas[0] = -0.99;

	grid->epsilons = Calloc(double, grid->Ne);
	grid->epsilons[0] = 1e-6;

	grid->weight_idxs = Calloc(double, grid->Nw);
	grid->weight_idxs[0] = 1;

	// start test code //
	struct GenQueue *q = gensvm_init_queue();
	gensvm_fill_queue(grid, q, train_data, test_data);

	mu_assert(q->N == 6, "Incorrect number of queue elements");

	int i;
	for (i=0; i<q->N; i++) {
		mu_assert(q->tasks[i]->ID == i, "Incorrect ID");
		mu_assert(q->tasks[i]->folds == 10, "Incorrect folds");
		mu_assert(q->tasks[i]->kerneltype == K_LINEAR,
				"Incorrect kernel type");
		mu_assert(q->tasks[i]->train_data == train_data,
				"Incorrect train_data");
		mu_assert(q->tasks[i]->test_data == test_data,
				"Incorrect test data")
	}

	// test p value
	mu_assert(q->tasks[0]->p == 1.0, "Incorrect p at task 0");
	mu_assert(q->tasks[1]->p == 1.5, "Incorrect p at task 1");
	mu_assert(q->tasks[2]->p == 2.0, "Incorrect p at task 2");
	mu_assert(q->tasks[3]->p == 1.0, "Incorrect p at task 3");
	mu_assert(q->tasks[4]->p == 1.5, "Incorrect p at task 4");
	mu_assert(q->tasks[5]->p == 2.0, "Incorrect p at task 5");

	// test lambda value
	mu_assert(q->tasks[0]->lambda == 1.0, "Incorrect lambda at task 0");
	mu_assert(q->tasks[1]->lambda == 1.0, "Incorrect lambda at task 1");
	mu_assert(q->tasks[2]->lambda == 1.0, "Incorrect lambda at task 2");
	mu_assert(q->tasks[3]->lambda == 5.0, "Incorrect lambda at task 3");
	mu_assert(q->tasks[4]->lambda == 5.0, "Incorrect lambda at task 4");
	mu_assert(q->tasks[5]->lambda == 5.0, "Incorrect lambda at task 5");

	// test kappa value
	mu_assert(q->tasks[0]->kappa == -0.99, "Incorrect kappa at task 0");
	mu_assert(q->tasks[1]->kappa == -0.99, "Incorrect kappa at task 1");
	mu_assert(q->tasks[2]->kappa == -0.99, "Incorrect kappa at task 2");
	mu_assert(q->tasks[3]->kappa == -0.99, "Incorrect kappa at task 3");
	mu_assert(q->tasks[4]->kappa == -0.99, "Incorrect kappa at task 4");
	mu_assert(q->tasks[5]->kappa == -0.99, "Incorrect kappa at task 5");

	// test epsilon value
	mu_assert(q->tasks[0]->epsilon == 1e-6, "Incorrect epsilon at task 0");
	mu_assert(q->tasks[1]->epsilon == 1e-6, "Incorrect epsilon at task 1");
	mu_assert(q->tasks[2]->epsilon == 1e-6, "Incorrect epsilon at task 2");
	mu_assert(q->tasks[3]->epsilon == 1e-6, "Incorrect epsilon at task 3");
	mu_assert(q->tasks[4]->epsilon == 1e-6, "Incorrect epsilon at task 4");
	mu_assert(q->tasks[5]->epsilon == 1e-6, "Incorrect epsilon at task 5");

	gensvm_free_queue(q);
	// end test code //

	gensvm_free_data(train_data);
	gensvm_free_data(test_data);
	gensvm_free_grid(grid);

	return NULL;
}

char *test_fill_queue_kernel()
{
	struct GenData *train_data = gensvm_init_data();
	struct GenData *test_data = gensvm_init_data();
	struct GenGrid *grid = gensvm_init_grid();

	grid->kerneltype = K_POLY;
	grid->Np = 3;
	grid->Nl = 2;
	grid->Nk = 1;
	grid->Ne = 1;
	grid->Nw = 1;
	grid->Ng = 2;
	grid->Nc = 3;
	grid->Nd = 2;

	grid->ps = Calloc(double, grid->Np);
	grid->ps[0] = 1.0;
	grid->ps[1] = 1.5;
	grid->ps[2] = 2.0;

	grid->lambdas = Calloc(double, grid->Nl);
	grid->lambdas[0] = 1.0;
	grid->lambdas[1] = 5.0;

	grid->kappas = Calloc(double, grid->Nk);
	grid->kappas[0] = -0.99;

	grid->epsilons = Calloc(double, grid->Ne);
	grid->epsilons[0] = 1e-6;

	grid->weight_idxs = Calloc(double, grid->Nw);
	grid->weight_idxs[0] = 1;

	grid->gammas = Calloc(double, grid->Ng);
	grid->gammas[0] = 0.5;
	grid->gammas[1] = 1.5;

	grid->coefs = Calloc(double, grid->Nc);
	grid->coefs[0] = 7.0;
	grid->coefs[1] = 11.0;
	grid->coefs[2] = 13.0;

	grid->degrees = Calloc(double, grid->Nd);
	grid->degrees[0] = 3.0;
	grid->degrees[1] = 5.0;

	// start test code //
	struct GenQueue *q = gensvm_init_queue();
	gensvm_fill_queue(grid, q, train_data, test_data);

	mu_assert(q->N == 72, "Incorrect number of queue elements");

	int i;
	for (i=0; i<q->N; i++) {
		mu_assert(q->tasks[i]->ID == i, "Incorrect ID");
		mu_assert(q->tasks[i]->folds == 10, "Incorrect folds");
		mu_assert(q->tasks[i]->kerneltype == K_POLY,
				"Incorrect kernel type");
		mu_assert(q->tasks[i]->train_data == train_data,
				"Incorrect train_data");
		mu_assert(q->tasks[i]->test_data == test_data,
				"Incorrect test data");
		mu_assert(q->tasks[i]->weight_idx == 1,
				"Incorrect weight idx");
		mu_assert(q->tasks[i]->epsilon == 1e-6,
				"Incorrect epsilon");
		mu_assert(q->tasks[i]->kappa == -0.99,
				"Incorrect kappa");
	}

	mu_assert(q->tasks[0]->p == 1.000000, "Incorrect p at task 0");
	mu_assert(q->tasks[1]->p == 1.500000, "Incorrect p at task 1");
	mu_assert(q->tasks[2]->p == 2.000000, "Incorrect p at task 2");
	mu_assert(q->tasks[3]->p == 1.000000, "Incorrect p at task 3");
	mu_assert(q->tasks[4]->p == 1.500000, "Incorrect p at task 4");
	mu_assert(q->tasks[5]->p == 2.000000, "Incorrect p at task 5");
	mu_assert(q->tasks[6]->p == 1.000000, "Incorrect p at task 6");
	mu_assert(q->tasks[7]->p == 1.500000, "Incorrect p at task 7");
	mu_assert(q->tasks[8]->p == 2.000000, "Incorrect p at task 8");
	mu_assert(q->tasks[9]->p == 1.000000, "Incorrect p at task 9");
	mu_assert(q->tasks[10]->p == 1.500000, "Incorrect p at task 10");
	mu_assert(q->tasks[11]->p == 2.000000, "Incorrect p at task 11");
	mu_assert(q->tasks[12]->p == 1.000000, "Incorrect p at task 12");
	mu_assert(q->tasks[13]->p == 1.500000, "Incorrect p at task 13");
	mu_assert(q->tasks[14]->p == 2.000000, "Incorrect p at task 14");
	mu_assert(q->tasks[15]->p == 1.000000, "Incorrect p at task 15");
	mu_assert(q->tasks[16]->p == 1.500000, "Incorrect p at task 16");
	mu_assert(q->tasks[17]->p == 2.000000, "Incorrect p at task 17");
	mu_assert(q->tasks[18]->p == 1.000000, "Incorrect p at task 18");
	mu_assert(q->tasks[19]->p == 1.500000, "Incorrect p at task 19");
	mu_assert(q->tasks[20]->p == 2.000000, "Incorrect p at task 20");
	mu_assert(q->tasks[21]->p == 1.000000, "Incorrect p at task 21");
	mu_assert(q->tasks[22]->p == 1.500000, "Incorrect p at task 22");
	mu_assert(q->tasks[23]->p == 2.000000, "Incorrect p at task 23");
	mu_assert(q->tasks[24]->p == 1.000000, "Incorrect p at task 24");
	mu_assert(q->tasks[25]->p == 1.500000, "Incorrect p at task 25");
	mu_assert(q->tasks[26]->p == 2.000000, "Incorrect p at task 26");
	mu_assert(q->tasks[27]->p == 1.000000, "Incorrect p at task 27");
	mu_assert(q->tasks[28]->p == 1.500000, "Incorrect p at task 28");
	mu_assert(q->tasks[29]->p == 2.000000, "Incorrect p at task 29");
	mu_assert(q->tasks[30]->p == 1.000000, "Incorrect p at task 30");
	mu_assert(q->tasks[31]->p == 1.500000, "Incorrect p at task 31");
	mu_assert(q->tasks[32]->p == 2.000000, "Incorrect p at task 32");
	mu_assert(q->tasks[33]->p == 1.000000, "Incorrect p at task 33");
	mu_assert(q->tasks[34]->p == 1.500000, "Incorrect p at task 34");
	mu_assert(q->tasks[35]->p == 2.000000, "Incorrect p at task 35");
	mu_assert(q->tasks[36]->p == 1.000000, "Incorrect p at task 36");
	mu_assert(q->tasks[37]->p == 1.500000, "Incorrect p at task 37");
	mu_assert(q->tasks[38]->p == 2.000000, "Incorrect p at task 38");
	mu_assert(q->tasks[39]->p == 1.000000, "Incorrect p at task 39");
	mu_assert(q->tasks[40]->p == 1.500000, "Incorrect p at task 40");
	mu_assert(q->tasks[41]->p == 2.000000, "Incorrect p at task 41");
	mu_assert(q->tasks[42]->p == 1.000000, "Incorrect p at task 42");
	mu_assert(q->tasks[43]->p == 1.500000, "Incorrect p at task 43");
	mu_assert(q->tasks[44]->p == 2.000000, "Incorrect p at task 44");
	mu_assert(q->tasks[45]->p == 1.000000, "Incorrect p at task 45");
	mu_assert(q->tasks[46]->p == 1.500000, "Incorrect p at task 46");
	mu_assert(q->tasks[47]->p == 2.000000, "Incorrect p at task 47");
	mu_assert(q->tasks[48]->p == 1.000000, "Incorrect p at task 48");
	mu_assert(q->tasks[49]->p == 1.500000, "Incorrect p at task 49");
	mu_assert(q->tasks[50]->p == 2.000000, "Incorrect p at task 50");
	mu_assert(q->tasks[51]->p == 1.000000, "Incorrect p at task 51");
	mu_assert(q->tasks[52]->p == 1.500000, "Incorrect p at task 52");
	mu_assert(q->tasks[53]->p == 2.000000, "Incorrect p at task 53");
	mu_assert(q->tasks[54]->p == 1.000000, "Incorrect p at task 54");
	mu_assert(q->tasks[55]->p == 1.500000, "Incorrect p at task 55");
	mu_assert(q->tasks[56]->p == 2.000000, "Incorrect p at task 56");
	mu_assert(q->tasks[57]->p == 1.000000, "Incorrect p at task 57");
	mu_assert(q->tasks[58]->p == 1.500000, "Incorrect p at task 58");
	mu_assert(q->tasks[59]->p == 2.000000, "Incorrect p at task 59");
	mu_assert(q->tasks[60]->p == 1.000000, "Incorrect p at task 60");
	mu_assert(q->tasks[61]->p == 1.500000, "Incorrect p at task 61");
	mu_assert(q->tasks[62]->p == 2.000000, "Incorrect p at task 62");
	mu_assert(q->tasks[63]->p == 1.000000, "Incorrect p at task 63");
	mu_assert(q->tasks[64]->p == 1.500000, "Incorrect p at task 64");
	mu_assert(q->tasks[65]->p == 2.000000, "Incorrect p at task 65");
	mu_assert(q->tasks[66]->p == 1.000000, "Incorrect p at task 66");
	mu_assert(q->tasks[67]->p == 1.500000, "Incorrect p at task 67");
	mu_assert(q->tasks[68]->p == 2.000000, "Incorrect p at task 68");
	mu_assert(q->tasks[69]->p == 1.000000, "Incorrect p at task 69");
	mu_assert(q->tasks[70]->p == 1.500000, "Incorrect p at task 70");
	mu_assert(q->tasks[71]->p == 2.000000, "Incorrect p at task 71");

	mu_assert(q->tasks[0]->lambda == 1.000000,
			"Incorrect lambda at task 0");
	mu_assert(q->tasks[1]->lambda == 1.000000,
			"Incorrect lambda at task 1");
	mu_assert(q->tasks[2]->lambda == 1.000000,
			"Incorrect lambda at task 2");
	mu_assert(q->tasks[3]->lambda == 5.000000,
			"Incorrect lambda at task 3");
	mu_assert(q->tasks[4]->lambda == 5.000000,
			"Incorrect lambda at task 4");
	mu_assert(q->tasks[5]->lambda == 5.000000,
			"Incorrect lambda at task 5");
	mu_assert(q->tasks[6]->lambda == 1.000000,
			"Incorrect lambda at task 6");
	mu_assert(q->tasks[7]->lambda == 1.000000,
			"Incorrect lambda at task 7");
	mu_assert(q->tasks[8]->lambda == 1.000000,
			"Incorrect lambda at task 8");
	mu_assert(q->tasks[9]->lambda == 5.000000,
			"Incorrect lambda at task 9");
	mu_assert(q->tasks[10]->lambda == 5.000000,
			"Incorrect lambda at task 10");
	mu_assert(q->tasks[11]->lambda == 5.000000,
			"Incorrect lambda at task 11");
	mu_assert(q->tasks[12]->lambda == 1.000000,
			"Incorrect lambda at task 12");
	mu_assert(q->tasks[13]->lambda == 1.000000,
			"Incorrect lambda at task 13");
	mu_assert(q->tasks[14]->lambda == 1.000000,
			"Incorrect lambda at task 14");
	mu_assert(q->tasks[15]->lambda == 5.000000,
			"Incorrect lambda at task 15");
	mu_assert(q->tasks[16]->lambda == 5.000000,
			"Incorrect lambda at task 16");
	mu_assert(q->tasks[17]->lambda == 5.000000,
			"Incorrect lambda at task 17");
	mu_assert(q->tasks[18]->lambda == 1.000000,
			"Incorrect lambda at task 18");
	mu_assert(q->tasks[19]->lambda == 1.000000,
			"Incorrect lambda at task 19");
	mu_assert(q->tasks[20]->lambda == 1.000000,
			"Incorrect lambda at task 20");
	mu_assert(q->tasks[21]->lambda == 5.000000,
			"Incorrect lambda at task 21");
	mu_assert(q->tasks[22]->lambda == 5.000000,
			"Incorrect lambda at task 22");
	mu_assert(q->tasks[23]->lambda == 5.000000,
			"Incorrect lambda at task 23");
	mu_assert(q->tasks[24]->lambda == 1.000000,
			"Incorrect lambda at task 24");
	mu_assert(q->tasks[25]->lambda == 1.000000,
			"Incorrect lambda at task 25");
	mu_assert(q->tasks[26]->lambda == 1.000000,
			"Incorrect lambda at task 26");
	mu_assert(q->tasks[27]->lambda == 5.000000,
			"Incorrect lambda at task 27");
	mu_assert(q->tasks[28]->lambda == 5.000000,
			"Incorrect lambda at task 28");
	mu_assert(q->tasks[29]->lambda == 5.000000,
			"Incorrect lambda at task 29");
	mu_assert(q->tasks[30]->lambda == 1.000000,
			"Incorrect lambda at task 30");
	mu_assert(q->tasks[31]->lambda == 1.000000,
			"Incorrect lambda at task 31");
	mu_assert(q->tasks[32]->lambda == 1.000000,
			"Incorrect lambda at task 32");
	mu_assert(q->tasks[33]->lambda == 5.000000,
			"Incorrect lambda at task 33");
	mu_assert(q->tasks[34]->lambda == 5.000000,
			"Incorrect lambda at task 34");
	mu_assert(q->tasks[35]->lambda == 5.000000,
			"Incorrect lambda at task 35");
	mu_assert(q->tasks[36]->lambda == 1.000000,
			"Incorrect lambda at task 36");
	mu_assert(q->tasks[37]->lambda == 1.000000,
			"Incorrect lambda at task 37");
	mu_assert(q->tasks[38]->lambda == 1.000000,
			"Incorrect lambda at task 38");
	mu_assert(q->tasks[39]->lambda == 5.000000,
			"Incorrect lambda at task 39");
	mu_assert(q->tasks[40]->lambda == 5.000000,
			"Incorrect lambda at task 40");
	mu_assert(q->tasks[41]->lambda == 5.000000,
			"Incorrect lambda at task 41");
	mu_assert(q->tasks[42]->lambda == 1.000000,
			"Incorrect lambda at task 42");
	mu_assert(q->tasks[43]->lambda == 1.000000,
			"Incorrect lambda at task 43");
	mu_assert(q->tasks[44]->lambda == 1.000000,
			"Incorrect lambda at task 44");
	mu_assert(q->tasks[45]->lambda == 5.000000,
			"Incorrect lambda at task 45");
	mu_assert(q->tasks[46]->lambda == 5.000000,
			"Incorrect lambda at task 46");
	mu_assert(q->tasks[47]->lambda == 5.000000,
			"Incorrect lambda at task 47");
	mu_assert(q->tasks[48]->lambda == 1.000000,
			"Incorrect lambda at task 48");
	mu_assert(q->tasks[49]->lambda == 1.000000,
			"Incorrect lambda at task 49");
	mu_assert(q->tasks[50]->lambda == 1.000000,
			"Incorrect lambda at task 50");
	mu_assert(q->tasks[51]->lambda == 5.000000,
			"Incorrect lambda at task 51");
	mu_assert(q->tasks[52]->lambda == 5.000000,
			"Incorrect lambda at task 52");
	mu_assert(q->tasks[53]->lambda == 5.000000,
			"Incorrect lambda at task 53");
	mu_assert(q->tasks[54]->lambda == 1.000000,
			"Incorrect lambda at task 54");
	mu_assert(q->tasks[55]->lambda == 1.000000,
			"Incorrect lambda at task 55");
	mu_assert(q->tasks[56]->lambda == 1.000000,
			"Incorrect lambda at task 56");
	mu_assert(q->tasks[57]->lambda == 5.000000,
			"Incorrect lambda at task 57");
	mu_assert(q->tasks[58]->lambda == 5.000000,
			"Incorrect lambda at task 58");
	mu_assert(q->tasks[59]->lambda == 5.000000,
			"Incorrect lambda at task 59");
	mu_assert(q->tasks[60]->lambda == 1.000000,
			"Incorrect lambda at task 60");
	mu_assert(q->tasks[61]->lambda == 1.000000,
			"Incorrect lambda at task 61");
	mu_assert(q->tasks[62]->lambda == 1.000000,
			"Incorrect lambda at task 62");
	mu_assert(q->tasks[63]->lambda == 5.000000,
			"Incorrect lambda at task 63");
	mu_assert(q->tasks[64]->lambda == 5.000000,
			"Incorrect lambda at task 64");
	mu_assert(q->tasks[65]->lambda == 5.000000,
			"Incorrect lambda at task 65");
	mu_assert(q->tasks[66]->lambda == 1.000000,
			"Incorrect lambda at task 66");
	mu_assert(q->tasks[67]->lambda == 1.000000,
			"Incorrect lambda at task 67");
	mu_assert(q->tasks[68]->lambda == 1.000000,
			"Incorrect lambda at task 68");
	mu_assert(q->tasks[69]->lambda == 5.000000,
			"Incorrect lambda at task 69");
	mu_assert(q->tasks[70]->lambda == 5.000000,
			"Incorrect lambda at task 70");
	mu_assert(q->tasks[71]->lambda == 5.000000,
			"Incorrect lambda at task 71");

	mu_assert(q->tasks[0]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 0");
	mu_assert(q->tasks[1]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 1");
	mu_assert(q->tasks[2]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 2");
	mu_assert(q->tasks[3]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 3");
	mu_assert(q->tasks[4]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 4");
	mu_assert(q->tasks[5]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 5");
	mu_assert(q->tasks[6]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 6");
	mu_assert(q->tasks[7]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 7");
	mu_assert(q->tasks[8]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 8");
	mu_assert(q->tasks[9]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 9");
	mu_assert(q->tasks[10]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 10");
	mu_assert(q->tasks[11]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 11");
	mu_assert(q->tasks[12]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 12");
	mu_assert(q->tasks[13]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 13");
	mu_assert(q->tasks[14]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 14");
	mu_assert(q->tasks[15]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 15");
	mu_assert(q->tasks[16]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 16");
	mu_assert(q->tasks[17]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 17");
	mu_assert(q->tasks[18]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 18");
	mu_assert(q->tasks[19]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 19");
	mu_assert(q->tasks[20]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 20");
	mu_assert(q->tasks[21]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 21");
	mu_assert(q->tasks[22]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 22");
	mu_assert(q->tasks[23]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 23");
	mu_assert(q->tasks[24]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 24");
	mu_assert(q->tasks[25]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 25");
	mu_assert(q->tasks[26]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 26");
	mu_assert(q->tasks[27]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 27");
	mu_assert(q->tasks[28]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 28");
	mu_assert(q->tasks[29]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 29");
	mu_assert(q->tasks[30]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 30");
	mu_assert(q->tasks[31]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 31");
	mu_assert(q->tasks[32]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 32");
	mu_assert(q->tasks[33]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 33");
	mu_assert(q->tasks[34]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 34");
	mu_assert(q->tasks[35]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 35");
	mu_assert(q->tasks[36]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 36");
	mu_assert(q->tasks[37]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 37");
	mu_assert(q->tasks[38]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 38");
	mu_assert(q->tasks[39]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 39");
	mu_assert(q->tasks[40]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 40");
	mu_assert(q->tasks[41]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 41");
	mu_assert(q->tasks[42]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 42");
	mu_assert(q->tasks[43]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 43");
	mu_assert(q->tasks[44]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 44");
	mu_assert(q->tasks[45]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 45");
	mu_assert(q->tasks[46]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 46");
	mu_assert(q->tasks[47]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 47");
	mu_assert(q->tasks[48]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 48");
	mu_assert(q->tasks[49]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 49");
	mu_assert(q->tasks[50]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 50");
	mu_assert(q->tasks[51]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 51");
	mu_assert(q->tasks[52]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 52");
	mu_assert(q->tasks[53]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 53");
	mu_assert(q->tasks[54]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 54");
	mu_assert(q->tasks[55]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 55");
	mu_assert(q->tasks[56]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 56");
	mu_assert(q->tasks[57]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 57");
	mu_assert(q->tasks[58]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 58");
	mu_assert(q->tasks[59]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 59");
	mu_assert(q->tasks[60]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 60");
	mu_assert(q->tasks[61]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 61");
	mu_assert(q->tasks[62]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 62");
	mu_assert(q->tasks[63]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 63");
	mu_assert(q->tasks[64]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 64");
	mu_assert(q->tasks[65]->kernelparam[0] == 0.500000,
			"Incorrect kernelparam 0 at task 65");
	mu_assert(q->tasks[66]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 66");
	mu_assert(q->tasks[67]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 67");
	mu_assert(q->tasks[68]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 68");
	mu_assert(q->tasks[69]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 69");
	mu_assert(q->tasks[70]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 70");
	mu_assert(q->tasks[71]->kernelparam[0] == 1.500000,
			"Incorrect kernelparam 0 at task 71");

	mu_assert(q->tasks[0]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 0");
	mu_assert(q->tasks[1]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 1");
	mu_assert(q->tasks[2]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 2");
	mu_assert(q->tasks[3]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 3");
	mu_assert(q->tasks[4]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 4");
	mu_assert(q->tasks[5]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 5");
	mu_assert(q->tasks[6]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 6");
	mu_assert(q->tasks[7]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 7");
	mu_assert(q->tasks[8]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 8");
	mu_assert(q->tasks[9]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 9");
	mu_assert(q->tasks[10]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 10");
	mu_assert(q->tasks[11]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 11");
	mu_assert(q->tasks[12]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 12");
	mu_assert(q->tasks[13]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 13");
	mu_assert(q->tasks[14]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 14");
	mu_assert(q->tasks[15]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 15");
	mu_assert(q->tasks[16]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 16");
	mu_assert(q->tasks[17]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 17");
	mu_assert(q->tasks[18]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 18");
	mu_assert(q->tasks[19]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 19");
	mu_assert(q->tasks[20]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 20");
	mu_assert(q->tasks[21]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 21");
	mu_assert(q->tasks[22]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 22");
	mu_assert(q->tasks[23]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 23");
	mu_assert(q->tasks[24]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 24");
	mu_assert(q->tasks[25]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 25");
	mu_assert(q->tasks[26]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 26");
	mu_assert(q->tasks[27]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 27");
	mu_assert(q->tasks[28]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 28");
	mu_assert(q->tasks[29]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 29");
	mu_assert(q->tasks[30]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 30");
	mu_assert(q->tasks[31]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 31");
	mu_assert(q->tasks[32]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 32");
	mu_assert(q->tasks[33]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 33");
	mu_assert(q->tasks[34]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 34");
	mu_assert(q->tasks[35]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 35");
	mu_assert(q->tasks[36]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 36");
	mu_assert(q->tasks[37]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 37");
	mu_assert(q->tasks[38]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 38");
	mu_assert(q->tasks[39]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 39");
	mu_assert(q->tasks[40]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 40");
	mu_assert(q->tasks[41]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 41");
	mu_assert(q->tasks[42]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 42");
	mu_assert(q->tasks[43]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 43");
	mu_assert(q->tasks[44]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 44");
	mu_assert(q->tasks[45]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 45");
	mu_assert(q->tasks[46]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 46");
	mu_assert(q->tasks[47]->kernelparam[1] == 7.000000,
			"Incorrect kernelparam 1 at task 47");
	mu_assert(q->tasks[48]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 48");
	mu_assert(q->tasks[49]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 49");
	mu_assert(q->tasks[50]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 50");
	mu_assert(q->tasks[51]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 51");
	mu_assert(q->tasks[52]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 52");
	mu_assert(q->tasks[53]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 53");
	mu_assert(q->tasks[54]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 54");
	mu_assert(q->tasks[55]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 55");
	mu_assert(q->tasks[56]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 56");
	mu_assert(q->tasks[57]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 57");
	mu_assert(q->tasks[58]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 58");
	mu_assert(q->tasks[59]->kernelparam[1] == 11.000000,
			"Incorrect kernelparam 1 at task 59");
	mu_assert(q->tasks[60]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 60");
	mu_assert(q->tasks[61]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 61");
	mu_assert(q->tasks[62]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 62");
	mu_assert(q->tasks[63]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 63");
	mu_assert(q->tasks[64]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 64");
	mu_assert(q->tasks[65]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 65");
	mu_assert(q->tasks[66]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 66");
	mu_assert(q->tasks[67]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 67");
	mu_assert(q->tasks[68]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 68");
	mu_assert(q->tasks[69]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 69");
	mu_assert(q->tasks[70]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 70");
	mu_assert(q->tasks[71]->kernelparam[1] == 13.000000,
			"Incorrect kernelparam 1 at task 71");

	mu_assert(q->tasks[0]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 0");
	mu_assert(q->tasks[1]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 1");
	mu_assert(q->tasks[2]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 2");
	mu_assert(q->tasks[3]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 3");
	mu_assert(q->tasks[4]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 4");
	mu_assert(q->tasks[5]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 5");
	mu_assert(q->tasks[6]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 6");
	mu_assert(q->tasks[7]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 7");
	mu_assert(q->tasks[8]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 8");
	mu_assert(q->tasks[9]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 9");
	mu_assert(q->tasks[10]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 10");
	mu_assert(q->tasks[11]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 11");
	mu_assert(q->tasks[12]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 12");
	mu_assert(q->tasks[13]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 13");
	mu_assert(q->tasks[14]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 14");
	mu_assert(q->tasks[15]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 15");
	mu_assert(q->tasks[16]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 16");
	mu_assert(q->tasks[17]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 17");
	mu_assert(q->tasks[18]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 18");
	mu_assert(q->tasks[19]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 19");
	mu_assert(q->tasks[20]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 20");
	mu_assert(q->tasks[21]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 21");
	mu_assert(q->tasks[22]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 22");
	mu_assert(q->tasks[23]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 23");
	mu_assert(q->tasks[24]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 24");
	mu_assert(q->tasks[25]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 25");
	mu_assert(q->tasks[26]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 26");
	mu_assert(q->tasks[27]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 27");
	mu_assert(q->tasks[28]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 28");
	mu_assert(q->tasks[29]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 29");
	mu_assert(q->tasks[30]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 30");
	mu_assert(q->tasks[31]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 31");
	mu_assert(q->tasks[32]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 32");
	mu_assert(q->tasks[33]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 33");
	mu_assert(q->tasks[34]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 34");
	mu_assert(q->tasks[35]->kernelparam[2] == 3.000000,
			"Incorrect kernelparam 2 at task 35");
	mu_assert(q->tasks[36]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 36");
	mu_assert(q->tasks[37]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 37");
	mu_assert(q->tasks[38]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 38");
	mu_assert(q->tasks[39]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 39");
	mu_assert(q->tasks[40]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 40");
	mu_assert(q->tasks[41]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 41");
	mu_assert(q->tasks[42]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 42");
	mu_assert(q->tasks[43]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 43");
	mu_assert(q->tasks[44]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 44");
	mu_assert(q->tasks[45]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 45");
	mu_assert(q->tasks[46]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 46");
	mu_assert(q->tasks[47]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 47");
	mu_assert(q->tasks[48]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 48");
	mu_assert(q->tasks[49]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 49");
	mu_assert(q->tasks[50]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 50");
	mu_assert(q->tasks[51]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 51");
	mu_assert(q->tasks[52]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 52");
	mu_assert(q->tasks[53]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 53");
	mu_assert(q->tasks[54]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 54");
	mu_assert(q->tasks[55]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 55");
	mu_assert(q->tasks[56]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 56");
	mu_assert(q->tasks[57]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 57");
	mu_assert(q->tasks[58]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 58");
	mu_assert(q->tasks[59]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 59");
	mu_assert(q->tasks[60]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 60");
	mu_assert(q->tasks[61]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 61");
	mu_assert(q->tasks[62]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 62");
	mu_assert(q->tasks[63]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 63");
	mu_assert(q->tasks[64]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 64");
	mu_assert(q->tasks[65]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 65");
	mu_assert(q->tasks[66]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 66");
	mu_assert(q->tasks[67]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 67");
	mu_assert(q->tasks[68]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 68");
	mu_assert(q->tasks[69]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 69");
	mu_assert(q->tasks[70]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 70");
	mu_assert(q->tasks[71]->kernelparam[2] == 5.000000,
			"Incorrect kernelparam 2 at task 71");
	gensvm_free_queue(q);
	// end test code //

	gensvm_free_data(train_data);
	gensvm_free_data(test_data);
	gensvm_free_grid(grid);

	return NULL;
}


char *test_kernel_changed()
{
	struct GenTask *new = gensvm_init_task();
	struct GenTask *old = gensvm_init_task();

	// start test code //
	mu_assert(gensvm_kernel_changed(new, NULL) == true,
			"Incorrect kernel changed (1)");

	mu_assert(gensvm_kernel_changed(new, old) == false,
			"Incorrect kernel changed (2)");

	new->kerneltype = K_RBF;
	mu_assert(gensvm_kernel_changed(new, old) == true,
			"Incorrect kernel changed (3)");

	// rbf kernel
	old->kerneltype = K_RBF;
	old->kernelparam = Malloc(double, 1);
	old->kernelparam[0] = 1.0;
	new->kernelparam = Malloc(double, 1);
	new->kernelparam[0] = 1.0;
	mu_assert(gensvm_kernel_changed(new, old) == false,
			"Incorrect kernel changed (4)");
	new->kernelparam[0] = 2.0;
	mu_assert(gensvm_kernel_changed(new, old) == true,
			"Incorrect kernel changed (5)");

	free(old->kernelparam);
	free(new->kernelparam);

	// poly kernel
	old->kerneltype = K_POLY;
	new->kerneltype = K_POLY;
	old->kernelparam = Malloc(double, 3);
	old->kernelparam[0] = 1.0;
	old->kernelparam[1] = 2.0;
	old->kernelparam[2] = 3.0;

	new->kernelparam = Malloc(double, 3);
	new->kernelparam[0] = 1.0;
	new->kernelparam[1] = 2.0;
	new->kernelparam[2] = 3.0;
	mu_assert(gensvm_kernel_changed(new, old) == false,
			"Incorrect kernel changed (6)");
	new->kernelparam[2] = 5.0;
	mu_assert(gensvm_kernel_changed(new, old) == true,
			"Incorrect kernel changed (7)");

	free(old->kernelparam);
	free(new->kernelparam);

	// sigmoid kernel
	old->kerneltype = K_SIGMOID;
	new->kerneltype = K_SIGMOID;
	old->kernelparam = Malloc(double, 2);
	old->kernelparam[0] = 1.0;
	old->kernelparam[1] = 2.0;
	new->kernelparam = Malloc(double, 2);
	new->kernelparam[0] = 1.0;
	new->kernelparam[1] = 2.0;

	mu_assert(gensvm_kernel_changed(new, old) == false,
			"Incorrect kernel changed (8)");
	new->kernelparam[1] = 5.0;
	mu_assert(gensvm_kernel_changed(new, old) == true,
			"Incorrect kernel changed (9)");

	// end test code //
	gensvm_free_task(new);
	gensvm_free_task(old);

	return NULL;
}

char *test_kernel_folds()
{
	mu_test_missing();

	return NULL;
}

char *test_train_queue()
{
	mu_test_missing();

	return NULL;
}

char *test_gridsearch_progress_linear()
{
	FILE *fid = NULL;
	const char *filename = "./data/test_progress_string.txt";
	GENSVM_OUTPUT_FILE = fopen(filename, "w");

	struct GenTask *task = gensvm_init_task();
	task->ID = 0;

	// start test code //
	gensvm_gridsearch_progress(task, 10, 0.5, 0.123, 0.7);
	fclose(GENSVM_OUTPUT_FILE);

	char buffer[MAX_LINE_LENGTH];
	fid = fopen(filename, "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	const char *expected = ("(001/010)\teps = 1e-06\tw = 1\tk = 0.00\t"
				"l = 1.000000\tp = 1.00\t\t0.500% (0.123s)\t"
				"(best = 0.700%)\n");
	mu_assert(strcmp(buffer, expected) == 0, "Incorrect progress string");

	fclose(fid);
	// end test code //
	gensvm_free_task(task);

	return NULL;
}

char *test_gridsearch_progress_rbf()
{
	FILE *fid = NULL;
	const char *filename = "./data/test_progress_string.txt";
	GENSVM_OUTPUT_FILE = fopen(filename, "w");

	struct GenTask *task = gensvm_init_task();
	task->ID = 0;
	task->kerneltype = K_RBF;
	task->kernelparam = Malloc(double, 1);
	task->kernelparam[0] = 3.0;

	// start test code //
	gensvm_gridsearch_progress(task, 10, 0.5, 0.123, 0.7);
	fclose(GENSVM_OUTPUT_FILE);

	char buffer[MAX_LINE_LENGTH];
	fid = fopen(filename, "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	const char *expected = ("(001/010)\tg = 3.000\teps = 1e-06\tw = 1\t"
			"k = 0.00\tl = 1.000000\tp = 1.00\t\t0.500% (0.123s)\t"
			"(best = 0.700%)\n");
	mu_assert(strcmp(buffer, expected) == 0, "Incorrect progress string");

	fclose(fid);
	// end test code //
	gensvm_free_task(task);

	return NULL;
}

char *test_gridsearch_progress_poly()
{
	FILE *fid = NULL;
	const char *filename = "./data/test_progress_string.txt";
	GENSVM_OUTPUT_FILE = fopen(filename, "w");

	struct GenTask *task = gensvm_init_task();
	task->ID = 0;
	task->kerneltype = K_POLY;
	task->kernelparam = Malloc(double, 3);
	task->kernelparam[0] = 3.0;
	task->kernelparam[1] = 1.0;
	task->kernelparam[2] = 2.0;

	// start test code //
	gensvm_gridsearch_progress(task, 10, 0.5, 0.123, 0.7);
	fclose(GENSVM_OUTPUT_FILE);

	char buffer[MAX_LINE_LENGTH];
	fid = fopen(filename, "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	const char *expected = ("(001/010)\t"
			"d = 2.00\tc = 1.00\tg = 3.000\t"
			"eps = 1e-06\tw = 1\tk = 0.00\t"
			"l = 1.000000\tp = 1.00\t\t0.500% (0.123s)\t"
			"(best = 0.700%)\n");
	mu_assert(strcmp(buffer, expected) == 0, "Incorrect progress string");

	fclose(fid);
	// end test code //
	gensvm_free_task(task);

	return NULL;
}

char *test_gridsearch_progress_sigmoid()
{
	FILE *fid = NULL;
	const char *filename = "./data/test_progress_string.txt";
	GENSVM_OUTPUT_FILE = fopen(filename, "w");

	struct GenTask *task = gensvm_init_task();
	task->ID = 0;
	task->kerneltype = K_SIGMOID;
	task->kernelparam = Malloc(double, 2);
	task->kernelparam[0] = 3.0;
	task->kernelparam[1] = 1.0;

	// start test code //
	gensvm_gridsearch_progress(task, 10, 0.5, 0.123, 0.7);
	fclose(GENSVM_OUTPUT_FILE);

	char buffer[MAX_LINE_LENGTH];
	fid = fopen(filename, "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	const char *expected = ("(001/010)\t"
			"c = 1.00\tg = 3.000\t"
			"eps = 1e-06\tw = 1\tk = 0.00\t"
			"l = 1.000000\tp = 1.00\t\t0.500% (0.123s)\t"
			"(best = 0.700%)\n");
	mu_assert(strcmp(buffer, expected) == 0, "Incorrect progress string");

	fclose(fid);
	// end test code //
	gensvm_free_task(task);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_fill_queue_nokernel);
	mu_run_test(test_fill_queue_kernel);
	mu_run_test(test_kernel_changed);
	mu_run_test(test_kernel_folds);
	mu_run_test(test_train_queue);
	mu_run_test(test_gridsearch_progress_linear);
	mu_run_test(test_gridsearch_progress_rbf);
	mu_run_test(test_gridsearch_progress_poly);
	mu_run_test(test_gridsearch_progress_sigmoid);

	return NULL;
}

RUN_TESTS(all_tests);
