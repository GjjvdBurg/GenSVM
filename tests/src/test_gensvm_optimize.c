/**
 * @file test_gensvm_optimize.c
 * @author G.J.J. van den Burg
 * @date 2016-09-01
 * @brief Unit tests for gensvm_optimize.c functions
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
#include "gensvm_optimize.h"
#include "gensvm_init.h"

#include "gensvm_debug.h"

extern FILE *GENSVM_OUTPUT_FILE;
extern FILE *GENSVM_ERROR_FILE;

char *test_gensvm_optimize()
{
	struct GenModel *model = gensvm_init_model();
	struct GenModel *seed_model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	int n = 8,
	    m = 3,
	    K = 4;
	data->n = n;
	data->m = m;
	data->r = m;
	data->K = K;

	model->n = n;
	model->m = m;
	model->K = K;

	seed_model->n = n;
	seed_model->m = m;
	seed_model->K = K;

	data->Z = Malloc(double, n*(m+1));
	data->y = Malloc(long, n);

	matrix_set(data->Z, data->m+1, 0, 0, 1.0);
	matrix_set(data->Z, data->m+1, 0, 1, 0.8740239771176158);
	matrix_set(data->Z, data->m+1, 0, 2, 0.3231542341162253);
	matrix_set(data->Z, data->m+1, 0, 3, 0.2533980609669184);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0);
	matrix_set(data->Z, data->m+1, 1, 1, 0.3433368959379667);
	matrix_set(data->Z, data->m+1, 1, 2, 0.2945713387329698);
	matrix_set(data->Z, data->m+1, 1, 3, 0.3042498181639990);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0);
	matrix_set(data->Z, data->m+1, 2, 1, 0.6513609117457242);
	matrix_set(data->Z, data->m+1, 2, 2, 0.7738077314847138);
	matrix_set(data->Z, data->m+1, 2, 3, 0.4426344045213226);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0);
	matrix_set(data->Z, data->m+1, 3, 1, 0.7223733317092962);
	matrix_set(data->Z, data->m+1, 3, 2, 0.9718611208972370);
	matrix_set(data->Z, data->m+1, 3, 3, 0.0796059591969125);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0);
	matrix_set(data->Z, data->m+1, 4, 1, 0.3014806706103061);
	matrix_set(data->Z, data->m+1, 4, 2, 0.1728058294642182);
	matrix_set(data->Z, data->m+1, 4, 3, 0.0851401652628196);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0);
	matrix_set(data->Z, data->m+1, 5, 1, 0.5114600128301799);
	matrix_set(data->Z, data->m+1, 5, 2, 0.3319865781913825);
	matrix_set(data->Z, data->m+1, 5, 3, 0.3330906711041684);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0);
	matrix_set(data->Z, data->m+1, 6, 1, 0.5824718351045201);
	matrix_set(data->Z, data->m+1, 6, 2, 0.7224023004247955);
	matrix_set(data->Z, data->m+1, 6, 3, 0.0937250920308128);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0);
	matrix_set(data->Z, data->m+1, 7, 1, 0.8228264179835741);
	matrix_set(data->Z, data->m+1, 7, 2, 0.4580785175957617);
	matrix_set(data->Z, data->m+1, 7, 3, 0.7585636149680212);

	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;
	data->y[5] = 2;
	data->y[6] = 4;
	data->y[7] = 1;

	model->p = 1.2143;
	model->kappa = 0.90298;
	model->lambda = 0.00219038;
	model->epsilon = 1e-15;

	gensvm_allocate_model(model);
	gensvm_allocate_model(seed_model);
	matrix_set(seed_model->V, K-1, 0, 0, 0.3294151808829250);
	matrix_set(seed_model->V, K-1, 0, 1, 0.8400578887926284);
	matrix_set(seed_model->V, K-1, 0, 2, 0.9336268164013294);
	matrix_set(seed_model->V, K-1, 1, 0, 0.6047157463292797);
	matrix_set(seed_model->V, K-1, 1, 1, 0.1390735925868357);
	matrix_set(seed_model->V, K-1, 1, 2, 0.6579825380479839);
	matrix_set(seed_model->V, K-1, 2, 0, 0.7628723943431572);
	matrix_set(seed_model->V, K-1, 2, 1, 0.3505528063594583);
	matrix_set(seed_model->V, K-1, 2, 2, 0.1221488022463632);
	matrix_set(seed_model->V, K-1, 3, 0, 0.4561071643209315);
	matrix_set(seed_model->V, K-1, 3, 1, 0.0840834388268874);
	matrix_set(seed_model->V, K-1, 3, 2, 0.5312457860071739);

	gensvm_init_V(seed_model, model, data);
	gensvm_initialize_weights(data, model);

	model->rho[0] = 0.3294151808829250;
	model->rho[1] = 0.6047157463292797;
	model->rho[2] = 0.7628723943431572;
	model->rho[3] = 0.4561071643209315;
	model->rho[4] = 0.8400578887926284;
	model->rho[5] = 0.1390735925868357;
	model->rho[6] = 0.3505528063594583;
	model->rho[7] = 0.0840834388268874;

	// start test code //
	GENSVM_OUTPUT_FILE = stdout;
	gensvm_optimize(model, data);

	gensvm_print_matrix(model->V, model->m+1, model->K-1);

	// end test code //

	gensvm_free_data(data);
	gensvm_free_model(model);
	gensvm_free_model(seed_model);

	mu_test_missing();
	return NULL;
}

char *test_gensvm_get_loss_1()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	int n = 8,
	    m = 3,
	    K = 3;
	model->n = n;
	model->m = m;
	model->K = K;
	struct GenWork *work = gensvm_init_work(model);

	// initialize the data
	data->n = n;
	data->K = K;
	data->m = m;


	data->y = Calloc(long, data->n);
	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;
	data->y[5] = 3;
	data->y[6] = 1;
	data->y[7] = 2;

	data->Z = Calloc(double, (data->n)*(data->m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.6112542725178001);
	matrix_set(data->Z, data->m+1, 0, 2, -0.7672096202890778);
	matrix_set(data->Z, data->m+1, 0, 3, -0.2600867145849611);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.5881180210361963);
	matrix_set(data->Z, data->m+1, 1, 2, -0.5419496202623567);
	matrix_set(data->Z, data->m+1, 1, 3, 0.7079932865564023);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, -0.9411484777876639);
	matrix_set(data->Z, data->m+1, 2, 2, -0.0251648291772256);
	matrix_set(data->Z, data->m+1, 2, 3, 0.5335722872738475);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, -0.6506872332924795);
	matrix_set(data->Z, data->m+1, 3, 2, -0.6277901989029552);
	matrix_set(data->Z, data->m+1, 3, 3, -0.1196037902922388);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, -0.9955402476800429);
	matrix_set(data->Z, data->m+1, 4, 2, -0.9514564047869466);
	matrix_set(data->Z, data->m+1, 4, 3, -0.1093968234456487);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, 0.3277661334163890);
	matrix_set(data->Z, data->m+1, 5, 2, 0.8271472175263959);
	matrix_set(data->Z, data->m+1, 5, 3, 0.6938788574898458);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, -0.8459013990907077);
	matrix_set(data->Z, data->m+1, 6, 2, -0.2453035880572786);
	matrix_set(data->Z, data->m+1, 6, 3, 0.0078257345629504);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, -0.4629532094536982);
	matrix_set(data->Z, data->m+1, 7, 2, 0.2935215202707828);
	matrix_set(data->Z, data->m+1, 7, 3, 0.0540516162042732);

	// initialize the model
	model->weight_idx = 1;
	model->kappa = 0.5;
	model->p = 1.5;
	model->lambda = 0.123;

	gensvm_allocate_model(model);
	gensvm_initialize_weights(data, model);
	gensvm_simplex(model);
	gensvm_simplex_diff(model);

	matrix_set(model->V, model->K-1, 0, 0, 0.6019309459245683);
	matrix_set(model->V, model->K-1, 0, 1, 0.0063825200426701);
	matrix_set(model->V, model->K-1, 1, 0, -0.9130102529085783);
	matrix_set(model->V, model->K-1, 1, 1, -0.8230766493212237);
	matrix_set(model->V, model->K-1, 2, 0, 0.5727079522160434);
	matrix_set(model->V, model->K-1, 2, 1, 0.6466468145039965);
	matrix_set(model->V, model->K-1, 3, 0, -0.8065680884346328);
	matrix_set(model->V, model->K-1, 3, 1, 0.5912336906588613);

	// start test code //
	double loss = gensvm_get_loss(model, data, work);

	mu_assert(fabs(loss - 0.903071383013108) < 1e-14,
			"Incorrect value of the loss");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);
	gensvm_free_work(work);


	return NULL;
}

char *test_gensvm_get_loss_2()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	int n = 8,
	    m = 3,
	    K = 3;
	model->n = n;
	model->m = m;
	model->K = K;
	struct GenWork *work = gensvm_init_work(model);

	// initialize the data
	data->n = n;
	data->K = K;
	data->m = m;

	data->y = Calloc(long, data->n);
	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;
	data->y[5] = 3;
	data->y[6] = 1;
	data->y[7] = 2;

	data->Z = Calloc(double, (data->n)*(data->m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.6112542725178001);
	matrix_set(data->Z, data->m+1, 0, 2, -0.7672096202890778);
	matrix_set(data->Z, data->m+1, 0, 3, -0.2600867145849611);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.5881180210361963);
	matrix_set(data->Z, data->m+1, 1, 2, -0.5419496202623567);
	matrix_set(data->Z, data->m+1, 1, 3, 0.7079932865564023);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, -0.9411484777876639);
	matrix_set(data->Z, data->m+1, 2, 2, -0.0251648291772256);
	matrix_set(data->Z, data->m+1, 2, 3, 0.5335722872738475);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, -0.6506872332924795);
	matrix_set(data->Z, data->m+1, 3, 2, -0.6277901989029552);
	matrix_set(data->Z, data->m+1, 3, 3, -0.1196037902922388);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, -0.9955402476800429);
	matrix_set(data->Z, data->m+1, 4, 2, -0.9514564047869466);
	matrix_set(data->Z, data->m+1, 4, 3, -0.1093968234456487);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, 0.3277661334163890);
	matrix_set(data->Z, data->m+1, 5, 2, 0.8271472175263959);
	matrix_set(data->Z, data->m+1, 5, 3, 0.6938788574898458);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, -0.8459013990907077);
	matrix_set(data->Z, data->m+1, 6, 2, -0.2453035880572786);
	matrix_set(data->Z, data->m+1, 6, 3, 0.0078257345629504);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, -0.4629532094536982);
	matrix_set(data->Z, data->m+1, 7, 2, 0.2935215202707828);
	matrix_set(data->Z, data->m+1, 7, 3, 0.0540516162042732);

	// initialize the model
	model->weight_idx = 2;
	model->kappa = 0.5;
	model->p = 1.5;
	model->lambda = 0.123;

	gensvm_allocate_model(model);
	gensvm_initialize_weights(data, model);
	gensvm_simplex(model);
	gensvm_simplex_diff(model);

	matrix_set(model->V, model->K-1, 0, 0, 0.6019309459245683);
	matrix_set(model->V, model->K-1, 0, 1, 0.0063825200426701);
	matrix_set(model->V, model->K-1, 1, 0, -0.9130102529085783);
	matrix_set(model->V, model->K-1, 1, 1, -0.8230766493212237);
	matrix_set(model->V, model->K-1, 2, 0, 0.5727079522160434);
	matrix_set(model->V, model->K-1, 2, 1, 0.6466468145039965);
	matrix_set(model->V, model->K-1, 3, 0, -0.8065680884346328);
	matrix_set(model->V, model->K-1, 3, 1, 0.5912336906588613);

	// start test code //
	double loss = gensvm_get_loss(model, data, work);

	mu_assert(fabs(loss - 0.972847045993281) < 1e-14,
			"Incorrect value of the loss");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);
	gensvm_free_work(work);

	return NULL;
}

char *test_gensvm_calculate_errors()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	int n = 8,
	    m = 3,
	    K = 3;

	// initialize the data
	data->n = n;
	data->K = K;
	data->m = m;

	data->y = Calloc(long, data->n);
	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;
	data->y[5] = 3;
	data->y[6] = 1;
	data->y[7] = 2;

	data->Z = Calloc(double, (data->n)*(data->m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.6112542725178001);
	matrix_set(data->Z, data->m+1, 0, 2, -0.7672096202890778);
	matrix_set(data->Z, data->m+1, 0, 3, -0.2600867145849611);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.5881180210361963);
	matrix_set(data->Z, data->m+1, 1, 2, -0.5419496202623567);
	matrix_set(data->Z, data->m+1, 1, 3, 0.7079932865564023);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, -0.9411484777876639);
	matrix_set(data->Z, data->m+1, 2, 2, -0.0251648291772256);
	matrix_set(data->Z, data->m+1, 2, 3, 0.5335722872738475);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, -0.6506872332924795);
	matrix_set(data->Z, data->m+1, 3, 2, -0.6277901989029552);
	matrix_set(data->Z, data->m+1, 3, 3, -0.1196037902922388);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, -0.9955402476800429);
	matrix_set(data->Z, data->m+1, 4, 2, -0.9514564047869466);
	matrix_set(data->Z, data->m+1, 4, 3, -0.1093968234456487);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, 0.3277661334163890);
	matrix_set(data->Z, data->m+1, 5, 2, 0.8271472175263959);
	matrix_set(data->Z, data->m+1, 5, 3, 0.6938788574898458);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, -0.8459013990907077);
	matrix_set(data->Z, data->m+1, 6, 2, -0.2453035880572786);
	matrix_set(data->Z, data->m+1, 6, 3, 0.0078257345629504);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, -0.4629532094536982);
	matrix_set(data->Z, data->m+1, 7, 2, 0.2935215202707828);
	matrix_set(data->Z, data->m+1, 7, 3, 0.0540516162042732);

	// initialize the model
	model->n = n;
	model->m = m;
	model->K = K;
	gensvm_allocate_model(model);
	gensvm_simplex(model);
	gensvm_simplex_diff(model);

	matrix_set(model->V, model->K-1, 0, 0, 0.6019309459245683);
	matrix_set(model->V, model->K-1, 0, 1, 0.0063825200426701);
	matrix_set(model->V, model->K-1, 1, 0, -0.9130102529085783);
	matrix_set(model->V, model->K-1, 1, 1, -0.8230766493212237);
	matrix_set(model->V, model->K-1, 2, 0, 0.5727079522160434);
	matrix_set(model->V, model->K-1, 2, 1, 0.6466468145039965);
	matrix_set(model->V, model->K-1, 3, 0, -0.8065680884346328);
	matrix_set(model->V, model->K-1, 3, 1, 0.5912336906588613);

	// start test code //
	double *ZV = Calloc(double, (data->n)*(data->K-1));
	gensvm_calculate_errors(model, data, ZV);

	// test ZV values
	mu_assert(fabs(matrix_get(ZV, data->K-1, 0, 0) -
				-0.1857598783645273) < 1e-14,
			"Incorrect value of ZV at 0, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 0, 1) -
				-1.1466122836367203) < 1e-14,
			"Incorrect value of ZV at 0, 1");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 1, 0) -
				-0.8164504861888491) < 1e-14,
			"Incorrect value of ZV at 1, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 1, 1) -
				-0.4095442019090963) < 1e-14,
			"Incorrect value of ZV at 1, 1");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 2, 0) -
				1.0164346780798894) < 1e-14,
			"Incorrect value of ZV at 2, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 2, 1) -
				1.0802130116671276) < 1e-14,
			"Incorrect value of ZV at 2, 1");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 3, 0) -
				0.9329432226278516) < 1e-14,
			"Incorrect value of ZV at 3, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 3, 1) -
				0.0652756651284464) < 1e-14,
			"Incorrect value of ZV at 3, 1");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 4, 0) -
				1.0541987367985999) < 1e-14,
			"Incorrect value of ZV at 4, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 4, 1) -
				0.1458531104005502) < 1e-14,
			"Incorrect value of ZV at 4, 1");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 5, 0) -
				0.2167303509991526) < 1e-14,
			"Incorrect value of ZV at 5, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 5, 1) -
				0.6817225403124993) < 1e-14,
			"Incorrect value of ZV at 5, 1");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 6, 0) -
				1.2274482928895278) < 1e-14,
			"Incorrect value of ZV at 6, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 6, 1) -
				0.5486262633865150) < 1e-14,
			"Incorrect value of ZV at 6, 1");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 7, 0) -
				1.1491177728196642) < 1e-14,
			"Incorrect value of ZV at 7, 0");
	mu_assert(fabs(matrix_get(ZV, data->K-1, 7, 1) -
				0.6091903890783275) < 1e-14,
			"Incorrect value of ZV at 7, 1");

	// test Q values
	mu_assert(fabs(matrix_get(model->Q, data->K, 0, 0) -
				-0.1857598783645273) < 1e-14,
			"Incorrect value of Q at 0, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 0, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 0, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 0, 2) -
				0.9001154267384245) < 1e-14,
			"Incorrect value of Q at 0, 2");
	mu_assert(fabs(matrix_get(model->Q, data->K, 1, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 1, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 1, 1) -
				0.8164504861888491) < 1e-14,
			"Incorrect value of Q at 1, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 1, 2) -
				0.7629009259203254) < 1e-14,
			"Incorrect value of Q at 1, 2");
	mu_assert(fabs(matrix_get(model->Q, data->K, 2, 0) -
				1.4437092486421736) < 1e-14,
			"Incorrect value of Q at 2, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 2, 1) -
				0.4272745705622841) < 1e-14,
			"Incorrect value of Q at 2, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 2, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 2, 2");
	mu_assert(fabs(matrix_get(model->Q, data->K, 3, 0) -
				0.9329432226278516) < 1e-14,
			"Incorrect value of Q at 3, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 3, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 3, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 3, 2) -
				0.4099412270637651) < 1e-14,
			"Incorrect value of Q at 3, 2");
	mu_assert(fabs(matrix_get(model->Q, data->K, 4, 0) -
				0.6534118672271528) < 1e-14,
			"Incorrect value of Q at 4, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 4, 1) -
				-0.4007868695714470) < 1e-14,
			"Incorrect value of Q at 4, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 4, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 4, 2");
	mu_assert(fabs(matrix_get(model->Q, data->K, 5, 0) -
				0.6987542137426619) < 1e-14,
			"Incorrect value of Q at 5, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 5, 1) -
				0.4820238627435093) < 1e-14,
			"Incorrect value of Q at 5, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 5, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 5, 2");
	mu_assert(fabs(matrix_get(model->Q, data->K, 6, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 6, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 6, 1) -
				-1.2274482928895278) < 1e-14,
			"Incorrect value of Q at 6, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 6, 2) -
				-1.0888484277208184) < 1e-14,
			"Incorrect value of Q at 6, 2");
	mu_assert(fabs(matrix_get(model->Q, data->K, 7, 0) -
				1.1491177728196642) < 1e-14,
			"Incorrect value of Q at 7, 0");
	mu_assert(fabs(matrix_get(model->Q, data->K, 7, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of Q at 7, 1");
	mu_assert(fabs(matrix_get(model->Q, data->K, 7, 2) -
				0.0469845337266742) < 1e-14,
			"Incorrect value of Q at 7, 2");

	free(ZV);
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_calculate_huber()
{
	struct GenModel *model = gensvm_init_model();
	model->n = 5;
	model->m = 3;
	model->K = 3;
	model->kappa = 0.5;

	gensvm_allocate_model(model);

	matrix_set(model->Q, model->K, 0, 0, -0.3386242674244120);
	matrix_set(model->Q, model->K, 0, 1, 1.0828252163937386);
	matrix_set(model->Q, model->K, 0, 2, 0.9734009993634181);
	matrix_set(model->Q, model->K, 1, 0, 1.3744927461858576);
	matrix_set(model->Q, model->K, 1, 1, 1.8086820272988162);
	matrix_set(model->Q, model->K, 1, 2, 0.9587412628706828);
	matrix_set(model->Q, model->K, 2, 0, -0.0530412768492290);
	matrix_set(model->Q, model->K, 2, 1, 0.4026826962708268);
	matrix_set(model->Q, model->K, 2, 2, -1.9705914880746673);
	matrix_set(model->Q, model->K, 3, 0, -0.8749982403551375);
	matrix_set(model->Q, model->K, 3, 1, 1.3981525936474806);
	matrix_set(model->Q, model->K, 3, 2, 0.5845478158465323);
	matrix_set(model->Q, model->K, 4, 0, 0.9594104113136890);
	matrix_set(model->Q, model->K, 4, 1, -0.7058945833639207);
	matrix_set(model->Q, model->K, 4, 2, -1.8413342248272893);

	// start test code //
	gensvm_calculate_huber(model);

	mu_assert(fabs(matrix_get(model->H, model->K, 0, 0) -
				0.5973049764458478) < 1e-14,
			"Incorrect H at 0, 0");
	mu_assert(fabs(matrix_get(model->H, model->K, 0, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect H at 0, 1");
	mu_assert(fabs(matrix_get(model->H, model->K, 0, 2) -
				0.0002358356116216) < 1e-14,
			"Incorrect H at 0, 2");
	mu_assert(fabs(matrix_get(model->H, model->K, 1, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect H at 1, 0");
	mu_assert(fabs(matrix_get(model->H, model->K, 1, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect H at 1, 1");
	mu_assert(fabs(matrix_get(model->H, model->K, 1, 2) -
				0.0005674277965020) < 1e-14,
			"Incorrect H at 1, 2");
	mu_assert(fabs(matrix_get(model->H, model->K, 2, 0) -
				0.3696319769160848) < 1e-14,
			"Incorrect H at 2, 0");
	mu_assert(fabs(matrix_get(model->H, model->K, 2, 1) -
				0.1189293204447631) < 1e-14,
			"Incorrect H at 2, 1");
	mu_assert(fabs(matrix_get(model->H, model->K, 2, 2) -
				2.2205914880746676) < 1e-14,
			"Incorrect H at 2, 2");
	mu_assert(fabs(matrix_get(model->H, model->K, 3, 0) -
				1.1249982403551375) < 1e-14,
			"Incorrect H at 3, 0");
	mu_assert(fabs(matrix_get(model->H, model->K, 3, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect H at 3, 1");
	mu_assert(fabs(matrix_get(model->H, model->K, 3, 2) -
				0.0575335057726290) < 1e-14,
			"Incorrect H at 3, 2");
	mu_assert(fabs(matrix_get(model->H, model->K, 4, 0) -
				0.0005491715699080) < 1e-14,
			"Incorrect H at 4, 0");
	mu_assert(fabs(matrix_get(model->H, model->K, 4, 1) -
				0.9558945833639207) < 1e-14,
			"Incorrect H at 4, 1");
	mu_assert(fabs(matrix_get(model->H, model->K, 4, 2) -
				2.0913342248272890) < 1e-14,
			"Incorrect H at 4, 2");

	// start end code //

	gensvm_free_model(model);

	return NULL;
}

char *test_gensvm_step_doubling()
{
	struct GenModel *model = gensvm_init_model();
	model->n = 5;
	model->m = 3;
	model->K = 3;

	gensvm_allocate_model(model);

	// start test code //

	matrix_set(model->V, model->K-1, 0, 0, 0.4534886648299394);
	matrix_set(model->V, model->K-1, 0, 1, 0.0744278734246461);
	matrix_set(model->V, model->K-1, 1, 0, 0.3119251404109698);
	matrix_set(model->V, model->K-1, 1, 1, 0.4656439597720683);
	matrix_set(model->V, model->K-1, 2, 0, 0.2011718791962903);
	matrix_set(model->V, model->K-1, 2, 1, 0.6500799120482493);
	matrix_set(model->V, model->K-1, 3, 0, 0.4203721613186416);
	matrix_set(model->V, model->K-1, 3, 1, 0.3322561487796912);

	matrix_set(model->Vbar, model->K-1, 0, 0, 0.4716285304362394);
	matrix_set(model->Vbar, model->K-1, 0, 1, 0.9580963971307287);
	matrix_set(model->Vbar, model->K-1, 1, 0, 0.4665492786460124);
	matrix_set(model->Vbar, model->K-1, 1, 1, 0.7584128769293659);
	matrix_set(model->Vbar, model->K-1, 2, 0, 0.0694310200377222);
	matrix_set(model->Vbar, model->K-1, 2, 1, 0.9819492320913891);
	matrix_set(model->Vbar, model->K-1, 3, 0, 0.9308026068356582);
	matrix_set(model->Vbar, model->K-1, 3, 1, 0.1323286413371241);

	gensvm_step_doubling(model);

	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 0) -
				0.4353487992236394) < 1e-14,
			"Incorrect V at 0, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 1) -
				-0.8092406502814364) < 1e-14,
			"Incorrect V at 0, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 0) -
				0.1573010021759272) < 1e-14,
			"Incorrect V at 1, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 1) -
				0.1728750426147708) < 1e-14,
			"Incorrect V at 1, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 0) -
				0.3329127383548583) < 1e-14,
			"Incorrect V at 2, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 1) -
				0.3182105920051096) < 1e-14,
			"Incorrect V at 2, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 0) -
				-0.0900582841983750) < 1e-14,
			"Incorrect V at 3, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 1) -
				0.5321836562222582) < 1e-14,
			"Incorrect V at 3, 1");

	// end test code //

	gensvm_free_model(model);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();

	mu_run_test(test_gensvm_calculate_errors);
	mu_run_test(test_gensvm_get_loss_1);
	mu_run_test(test_gensvm_get_loss_2);
	mu_run_test(test_gensvm_calculate_huber);
	mu_run_test(test_gensvm_step_doubling);

	mu_run_test(test_gensvm_optimize);

	return NULL;
}

RUN_TESTS(all_tests);
