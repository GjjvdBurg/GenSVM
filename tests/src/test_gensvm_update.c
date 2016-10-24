/**
 * @file test_gensvm_update.c
 * @author G.J.J. van den Burg
 * @date 2016-09-01
 * @brief Unit tests for gensvm_update.c functions
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
#include "gensvm_simplex.h"

char *test_gensvm_calculate_omega()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	int n = 5,
	    m = 3,
	    K = 3;

	data->n = n;
	data->m = m;
	data->K = K;

	data->y = Calloc(long, n);
	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;

	model->n = n;
	model->m = m;
	model->K = K;
	model->p = 1.213;
	gensvm_allocate_model(model);

	matrix_set(model->H, model->K, 0, 0, 0.8465725800087526);
	matrix_set(model->H, model->K, 0, 1, 1.2876921677680249);
	matrix_set(model->H, model->K, 0, 2, 1.0338561593991831);
	matrix_set(model->H, model->K, 1, 0, 1.1891038526621391);
	matrix_set(model->H, model->K, 1, 1, 0.4034192031226095);
	matrix_set(model->H, model->K, 1, 2, 1.5298894170910078);
	matrix_set(model->H, model->K, 2, 0, 1.3505111116922732);
	matrix_set(model->H, model->K, 2, 1, 1.4336863304586636);
	matrix_set(model->H, model->K, 2, 2, 1.7847533480330757);
	matrix_set(model->H, model->K, 3, 0, 1.7712504341475415);
	matrix_set(model->H, model->K, 3, 1, 1.6905146737773038);
	matrix_set(model->H, model->K, 3, 2, 0.8189336598535132);
	matrix_set(model->H, model->K, 4, 0, 0.6164203008844277);
	matrix_set(model->H, model->K, 4, 1, 0.2456444285093894);
	matrix_set(model->H, model->K, 4, 2, 0.8184193969741095);

	// start test code //
	mu_assert(fabs(gensvm_calculate_omega(model, data, 0) -
				0.7394076262220608) < 1e-14,
			"Incorrect omega at 0");
	mu_assert(fabs(gensvm_calculate_omega(model, data, 1) -
				0.7294526264247443) < 1e-14,
			"Incorrect omega at 1");
	mu_assert(fabs(gensvm_calculate_omega(model, data, 2) -
				0.6802499471888741) < 1e-14,
			"Incorrect omega at 2");
	mu_assert(fabs(gensvm_calculate_omega(model, data, 3) -
				0.6886792032441273) < 1e-14,
			"Incorrect omega at 3");
	mu_assert(fabs(gensvm_calculate_omega(model, data, 4) -
				0.8695329737474283) < 1e-14,
			"Incorrect omega at 4");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_majorize_is_simple()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	int n = 5,
	    m = 3,
	    K = 3;

	data->n = n;
	data->m = m;
	data->K = K;

	data->y = Calloc(long, n);
	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;

	model->n = n;
	model->m = m;
	model->K = K;
	model->p = 1.213;
	gensvm_allocate_model(model);

	matrix_set(model->H, model->K, 0, 0, 0.8465725800087526);
	matrix_set(model->H, model->K, 0, 1, 1.2876921677680249);
	matrix_set(model->H, model->K, 0, 2, 1.0338561593991831);
	matrix_set(model->H, model->K, 1, 0, 1.1891038526621391);
	matrix_set(model->H, model->K, 1, 1, 0.4034192031226095);
	matrix_set(model->H, model->K, 1, 2, 0.0);
	matrix_set(model->H, model->K, 2, 0, 0.5);
	matrix_set(model->H, model->K, 2, 1, 0.0);
	matrix_set(model->H, model->K, 2, 2, 1.1);
	matrix_set(model->H, model->K, 3, 0, 0.0);
	matrix_set(model->H, model->K, 3, 1, 0.0);
	matrix_set(model->H, model->K, 3, 2, 0.8189336598535132);
	matrix_set(model->H, model->K, 4, 0, 0.6164203008844277);
	matrix_set(model->H, model->K, 4, 1, 0.2456444285093894);
	matrix_set(model->H, model->K, 4, 2, 0.8184193969741095);

	// start test code //
	mu_assert(gensvm_majorize_is_simple(model, data, 0) == false,
			"Incorrect simple at 0");
	mu_assert(gensvm_majorize_is_simple(model, data, 1) == true,
			"Incorrect simple at 1");
	mu_assert(gensvm_majorize_is_simple(model, data, 2) == true,
			"Incorrect simple at 2");
	mu_assert(gensvm_majorize_is_simple(model, data, 3) == true,
			"Incorrect simple at 3");
	mu_assert(gensvm_majorize_is_simple(model, data, 4) == false,
			"Incorrect simple at 4");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);
	return NULL;
}

char *test_gensvm_calculate_ab_non_simple()
{
	struct GenModel *model = gensvm_init_model();
	int n = 1,
	    m = 1,
	    K = 1;

	model->n = n;
	model->m = m;
	model->K = K;
	model->kappa = 0.5;

	gensvm_allocate_model(model);

	// start test code //
	double a, b_aq;

	// tests with p = 2
	model->p = 2.0;
	matrix_set(model->Q, K, 0, 0, -1.0);
	gensvm_calculate_ab_non_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 1.5) < 1e-14, "Incorrect value for a (1)");
	mu_assert(fabs(b_aq - 1.25) < 1e-14, "Incorrect value for b (1)");

	matrix_set(model->Q, K, 0, 0, 0.5);
	gensvm_calculate_ab_non_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 1.5) < 1e-14, "Incorrect value for a (2)");
	mu_assert(fabs(b_aq - 0.0277777777777778) < 1e-14,
			"Incorrect value for b (2)");

	matrix_set(model->Q, K, 0, 0, 2.0);
	gensvm_calculate_ab_non_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 1.5) < 1e-14, "Incorrect value for a (3)");
	mu_assert(fabs(b_aq - 0.0) < 1e-14, "Incorrect value for b (3)");

	// tests with p != 2 (namely, 1.5)
	// Note that here (p + kappa - 1)/(p - 2) = -2
	//
	// We distinguish: q <= (p + kappa - 1)/(p - 2)
	// 		   q in (p + kappa - 1)/(p - 2), -kappa]
	// 		   q in (-kappa, 1]
	// 		   q > 1
	model->p = 1.5;
	matrix_set(model->Q, K, 0, 0, -3.0);
	gensvm_calculate_ab_non_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 0.312018860376691) < 1e-14,
			"Incorrect value for a (4)");
	mu_assert(fabs(b_aq - 1.35208172829900) < 1e-14,
			"Incorrect value for b (4)");

	matrix_set(model->Q, K, 0, 0, -1.0);
	gensvm_calculate_ab_non_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 0.866025403784439) < 1e-14,
			"Incorrect value for a (5)");
	mu_assert(fabs(b_aq - 0.838525491562421) < 1e-14,
			"Incorrect value for b (5)");

	matrix_set(model->Q, K, 0, 0, 0.5);
	gensvm_calculate_ab_non_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 0.866025403784439) < 1e-14,
			"Incorrect value for a (6)");
	mu_assert(fabs(b_aq - 0.0721687836487032) < 1e-14,
			"Incorrect value for b (6)");

	matrix_set(model->Q, K, 0, 0, 2.0);
	gensvm_calculate_ab_non_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 0.245495126515491) < 1e-14,
			"Incorrect value for a (7)");
	mu_assert(fabs(b_aq - 0.0) < 1e-14,
			"Incorrect value for b (7)");
	// end test code //

	gensvm_free_model(model);

	return NULL;
}

char *test_gensvm_calculate_ab_simple()
{
	struct GenModel *model = gensvm_init_model();
	int n = 1,
	    m = 1,
	    K = 1;

	model->n = n;
	model->m = m;
	model->K = K;
	model->kappa = 0.5;

	gensvm_allocate_model(model);

	// start test code //
	double a, b_aq;

	matrix_set(model->Q, K, 0, 0, -1.5);
	gensvm_calculate_ab_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 0.142857142857143) < 1e-14,
			"Incorrect value for a (1)");
	mu_assert(fabs(b_aq - 0.5) < 1e-14,
			"Incorrect value for b (1)");

	matrix_set(model->Q, K, 0, 0, 0.75);
	gensvm_calculate_ab_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 0.333333333333333) < 1e-14,
			"Incorrect value for a (2)");
	mu_assert(fabs(b_aq - 0.0833333333333333) < 1e-14,
			"Incorrect value for b (2)");

	matrix_set(model->Q, K, 0, 0, 2.0);
	gensvm_calculate_ab_simple(model, 0, 0, &a, &b_aq);
	mu_assert(fabs(a - 0.142857142857143) < 1e-14,
			"Incorrect value for a (3)");
	mu_assert(fabs(b_aq - 0.0) < 1e-14,
			"Incorrect value for b (3)");

	// end test code //
	gensvm_free_model(model);

	return NULL;
}

char *test_gensvm_get_update()
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

	// initialize data
	data->n = n;
	data->m = m;
	data->K = K;

	data->y = Calloc(long, n);
	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;
	data->y[5] = 3;
	data->y[6] = 1;
	data->y[7] = 2;

	data->Z = Calloc(double, n*(m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.6437306339619082);
	matrix_set(data->Z, data->m+1, 0, 2, -0.3276778319121999);
	matrix_set(data->Z, data->m+1, 0, 3, 0.1564053473463392);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, -0.8683091763200105);
	matrix_set(data->Z, data->m+1, 1, 2, -0.6910830836015162);
	matrix_set(data->Z, data->m+1, 1, 3, -0.9675430665130734);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, -0.5024888699077029);
	matrix_set(data->Z, data->m+1, 2, 2, -0.9649738292750712);
	matrix_set(data->Z, data->m+1, 2, 3, 0.0776560791351473);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, 0.8206429991392579);
	matrix_set(data->Z, data->m+1, 3, 2, -0.7255681388968501);
	matrix_set(data->Z, data->m+1, 3, 3, -0.9475952272877165);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, 0.3426050950418613);
	matrix_set(data->Z, data->m+1, 4, 2, -0.5340602451864306);
	matrix_set(data->Z, data->m+1, 4, 3, -0.7159704241662815);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, -0.3077314049206620);
	matrix_set(data->Z, data->m+1, 5, 2, 0.1141288036288195);
	matrix_set(data->Z, data->m+1, 5, 3, -0.7060114827535847);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, 0.6301294373610109);
	matrix_set(data->Z, data->m+1, 6, 2, -0.9983027363627769);
	matrix_set(data->Z, data->m+1, 6, 3, -0.9365684178444004);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, -0.0665379368401439);
	matrix_set(data->Z, data->m+1, 7, 2, -0.1781385556871763);
	matrix_set(data->Z, data->m+1, 7, 3, -0.7292593770500276);

	// initialize model
	model->p = 1.1;
	model->lambda = 0.123;
	model->weight_idx = 1;
	model->kappa = 0.5;

	// initialize matrices
	gensvm_allocate_model(model);
	gensvm_initialize_weights(data, model);
	gensvm_simplex(model);
	gensvm_simplex_diff(model);

	// initialize V
	matrix_set(model->V, model->K-1, 0, 0, -0.7593642121025029);
	matrix_set(model->V, model->K-1, 0, 1, -0.5497320698504756);
	matrix_set(model->V, model->K-1, 1, 0, 0.2982680646268177);
	matrix_set(model->V, model->K-1, 1, 1, -0.2491408622891925);
	matrix_set(model->V, model->K-1, 2, 0, -0.3118572761092807);
	matrix_set(model->V, model->K-1, 2, 1, 0.5461219445756100);
	matrix_set(model->V, model->K-1, 3, 0, -0.3198994238626641);
	matrix_set(model->V, model->K-1, 3, 1, 0.7134997072555367);

	// start test code //

	// these need to be prepared for the update call
	gensvm_calculate_errors(model, data, work->ZV);
	gensvm_calculate_huber(model);

	// run the actual update call
	gensvm_get_update(model, data, work);

	// test values
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 0) -
				-0.1323791019594062) < 1e-14,
			"Incorrect value of model->V at 0, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 1) -
				-0.3598407983154332) < 1e-14,
			"Incorrect value of model->V at 0, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 0) -
				0.3532993103400935) < 1e-14,
			"Incorrect value of model->V at 1, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 1) -
				-0.4094572388475382) < 1e-14,
			"Incorrect value of model->V at 1, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 0) -
				0.1313169839871234) < 1e-14,
			"Incorrect value of model->V at 2, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 1) -
				0.2423439972728328) < 1e-14,
			"Incorrect value of model->V at 2, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 0) -
				0.0458431025455224) < 1e-14,
			"Incorrect value of model->V at 3, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 1) -
				0.4390030236354089) < 1e-14,
			"Incorrect value of model->V at 3, 1");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);
	gensvm_free_work(work);

	return NULL;
}

char *test_gensvm_get_update_sparse()
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

	// initialize data
	data->n = n;
	data->m = m;
	data->K = K;

	data->y = Calloc(long, n);
	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;
	data->y[5] = 3;
	data->y[6] = 1;
	data->y[7] = 2;

	data->Z = Calloc(double, n*(m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.6437306339619082);
	matrix_set(data->Z, data->m+1, 0, 2, -0.3276778319121999);
	matrix_set(data->Z, data->m+1, 0, 3, 0.1564053473463392);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, -0.8683091763200105);
	matrix_set(data->Z, data->m+1, 1, 2, -0.6910830836015162);
	matrix_set(data->Z, data->m+1, 1, 3, -0.9675430665130734);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, -0.5024888699077029);
	matrix_set(data->Z, data->m+1, 2, 2, -0.9649738292750712);
	matrix_set(data->Z, data->m+1, 2, 3, 0.0776560791351473);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, 0.8206429991392579);
	matrix_set(data->Z, data->m+1, 3, 2, -0.7255681388968501);
	matrix_set(data->Z, data->m+1, 3, 3, -0.9475952272877165);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, 0.3426050950418613);
	matrix_set(data->Z, data->m+1, 4, 2, -0.5340602451864306);
	matrix_set(data->Z, data->m+1, 4, 3, -0.7159704241662815);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, -0.3077314049206620);
	matrix_set(data->Z, data->m+1, 5, 2, 0.1141288036288195);
	matrix_set(data->Z, data->m+1, 5, 3, -0.7060114827535847);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, 0.6301294373610109);
	matrix_set(data->Z, data->m+1, 6, 2, -0.9983027363627769);
	matrix_set(data->Z, data->m+1, 6, 3, -0.9365684178444004);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, -0.0665379368401439);
	matrix_set(data->Z, data->m+1, 7, 2, -0.1781385556871763);
	matrix_set(data->Z, data->m+1, 7, 3, -0.7292593770500276);

	// convert Z to a sparse matrix to test the sparse functions
	data->spZ = gensvm_dense_to_sparse(data->Z, data->n, data->m+1);
	free(data->RAW);
	data->RAW = NULL;
	free(data->Z);
	data->Z = NULL;

	// initialize model
	model->p = 1.1;
	model->lambda = 0.123;
	model->weight_idx = 1;
	model->kappa = 0.5;

	// initialize matrices
	gensvm_allocate_model(model);
	gensvm_initialize_weights(data, model);
	gensvm_simplex(model);
	gensvm_simplex_diff(model);

	// initialize V
	matrix_set(model->V, model->K-1, 0, 0, -0.7593642121025029);
	matrix_set(model->V, model->K-1, 0, 1, -0.5497320698504756);
	matrix_set(model->V, model->K-1, 1, 0, 0.2982680646268177);
	matrix_set(model->V, model->K-1, 1, 1, -0.2491408622891925);
	matrix_set(model->V, model->K-1, 2, 0, -0.3118572761092807);
	matrix_set(model->V, model->K-1, 2, 1, 0.5461219445756100);
	matrix_set(model->V, model->K-1, 3, 0, -0.3198994238626641);
	matrix_set(model->V, model->K-1, 3, 1, 0.7134997072555367);

	// start test code //

	// these need to be prepared for the update call
	gensvm_calculate_errors(model, data, work->ZV);
	gensvm_calculate_huber(model);

	// run the actual update call
	gensvm_get_update(model, data, work);

	// test values
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 0) -
				-0.1323791019594062) < 1e-14,
			"Incorrect value of model->V at 0, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 1) -
				-0.3598407983154332) < 1e-14,
			"Incorrect value of model->V at 0, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 0) -
				0.3532993103400935) < 1e-14,
			"Incorrect value of model->V at 1, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 1) -
				-0.4094572388475382) < 1e-14,
			"Incorrect value of model->V at 1, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 0) -
				0.1313169839871234) < 1e-14,
			"Incorrect value of model->V at 2, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 1) -
				0.2423439972728328) < 1e-14,
			"Incorrect value of model->V at 2, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 0) -
				0.0458431025455224) < 1e-14,
			"Incorrect value of model->V at 3, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 1) -
				0.4390030236354089) < 1e-14,
			"Incorrect value of model->V at 3, 1");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);
	gensvm_free_work(work);

	return NULL;
}


char *test_dposv()
{
	int n = 6,
	    m = 5;

	// start test code //
	double *A = Calloc(double, n*n);
	double *B = Calloc(double, n*m);

	// We're only storing the upper triangular part of the symmetric
	// matrix A.
	matrix_set(A, n, 0, 0, 6.2522023496540386);
	matrix_set(A, n, 0, 1, 1.2760969977888754);
	matrix_set(A, n, 0, 2, 1.1267774552193974);
	matrix_set(A, n, 0, 3, 0.8384267227932789);
	matrix_set(A, n, 0, 4, 0.9588857509656767);
	matrix_set(A, n, 0, 5, 0.7965747978871199);
	matrix_set(A, n, 1, 1, 6.7539376310748924);
	matrix_set(A, n, 1, 2, 0.5908599276261999);
	matrix_set(A, n, 1, 3, 0.9987368128192129);
	matrix_set(A, n, 1, 4, 1.1142949385131484);
	matrix_set(A, n, 1, 5, 1.4150107613377123);
	matrix_set(A, n, 2, 2, 7.3361678639533139);
	matrix_set(A, n, 2, 3, 1.5596679563906113);
	matrix_set(A, n, 2, 4, 0.8162441257417704);
	matrix_set(A, n, 2, 5, 0.8701893160678078);
	matrix_set(A, n, 3, 3, 6.8330423955320834);
	matrix_set(A, n, 3, 4, 1.6081779105091201);
	matrix_set(A, n, 3, 5, 1.0498769837396527);
	matrix_set(A, n, 4, 4, 7.6810607313742949);
	matrix_set(A, n, 4, 5, 1.1352511350739003);
	matrix_set(A, n, 5, 5, 7.0573435553104567);

	// this is the matrix B (n x m), stored in COLUMN-MAJOR ORDER
	B[0] = 0.5759809004700531;
	B[1] = 0.4643751885289473;
	B[2] = 0.1914807543974765;
	B[3] = 0.2875503245961965;
	B[4] = 0.0493123646253395;
	B[5] = 0.4333053066976881;
	B[6] = 0.4738306027724854;
	B[7] = 0.2460182087225041;
	B[8] = 0.1620492662433550;
	B[9] = 0.9596380576403235;
	B[10] = 0.7244837218691289;
	B[11] = 0.9437116578537014;
	B[12] = 0.3320986772155025;
	B[13] = 0.4717424581951766;
	B[14] = 0.9206089360217588;
	B[15] = 0.7059004575000609;
	B[16] = 0.1696670763906902;
	B[17] = 0.4896586269167711;
	B[18] = 0.9539497766794410;
	B[19] = 0.2269749103273779;
	B[20] = 0.8832156948007016;
	B[21] = 0.4682217970327739;
	B[22] = 0.5293439096127632;
	B[23] = 0.8699136677253214;
	B[24] = 0.1622687366790325;
	B[25] = 0.4511598310105013;
	B[26] = 0.5587302139109592;
	B[27] = 0.7456952498557438;
	B[28] = 0.5923112589693547;
	B[29] = 0.2243481938151050;

	// note the 'L' here denotes the lower triangular part of A. Above we
	// stored the upper triangular part of A in row-major order, so that's
	// the lower triangular part in column-major order, which Lapack uses.
	int status = dposv('L', n, m, A, n, B, n);
	mu_assert(status == 0, "dposv didn't return status success");

	// Since B now contains the solution in Column-Major order, we have to
	// check it in that order.

	mu_assert(fabs(B[0] - 0.0770250502756885) < 1e-14,
			"Incorrect value of B at 0");
	mu_assert(fabs(B[1] - 0.0449013611583528) < 1e-14,
			"Incorrect value of B at 1");
	mu_assert(fabs(B[2] - 0.0028421256926631) < 1e-14,
			"Incorrect value of B at 2");
	mu_assert(fabs(B[3] - 0.0238082780914757) < 1e-14,
			"Incorrect value of B at 3");
	mu_assert(fabs(B[4] - -0.0213884392480803) < 1e-14,
			"Incorrect value of B at 4");
	mu_assert(fabs(B[5] - 0.0432493445363141) < 1e-14,
			"Incorrect value of B at 5");
	mu_assert(fabs(B[6] - 0.0469188408250497) < 1e-14,
			"Incorrect value of B at 6");
	mu_assert(fabs(B[7] - -0.0188676544565197) < 1e-14,
			"Incorrect value of B at 7");
	mu_assert(fabs(B[8] - -0.0268693544126544) < 1e-14,
			"Incorrect value of B at 8");
	mu_assert(fabs(B[9] - 0.1139942447258460) < 1e-14,
			"Incorrect value of B at 9");
	mu_assert(fabs(B[10] - 0.0539483576192093) < 1e-14,
			"Incorrect value of B at 10");
	mu_assert(fabs(B[11] - 0.1098843745987866) < 1e-14,
			"Incorrect value of B at 11");
	mu_assert(fabs(B[12] - 0.0142822905211067) < 1e-14,
			"Incorrect value of B at 12");
	mu_assert(fabs(B[13] - 0.0425078586146396) < 1e-14,
			"Incorrect value of B at 13");
	mu_assert(fabs(B[14] - 0.1022650353097420) < 1e-14,
			"Incorrect value of B at 14");
	mu_assert(fabs(B[15] - 0.0700476338859921) < 1e-14,
			"Incorrect value of B at 15");
	mu_assert(fabs(B[16] - -0.0171546096353451) < 1e-14,
			"Incorrect value of B at 16");
	mu_assert(fabs(B[17] - 0.0389772844468578) < 1e-14,
			"Incorrect value of B at 17");
	mu_assert(fabs(B[18] - 0.1231757430810565) < 1e-14,
			"Incorrect value of B at 18");
	mu_assert(fabs(B[19] - -0.0246954324681607) < 1e-14,
			"Incorrect value of B at 19");
	mu_assert(fabs(B[20] - 0.0853098528328778) < 1e-14,
			"Incorrect value of B at 20");
	mu_assert(fabs(B[21] - 0.0155226252622933) < 1e-14,
			"Incorrect value of B at 21");
	mu_assert(fabs(B[22] - 0.0305321945431931) < 1e-14,
			"Incorrect value of B at 22");
	mu_assert(fabs(B[23] - 0.0965724559730953) < 1e-14,
			"Incorrect value of B at 23");
	mu_assert(fabs(B[24] - -0.0106596940426243) < 1e-14,
			"Incorrect value of B at 24");
	mu_assert(fabs(B[25] - 0.0446093337039209) < 1e-14,
			"Incorrect value of B at 25");
	mu_assert(fabs(B[26] - 0.0517721408799503) < 1e-14,
			"Incorrect value of B at 26");
	mu_assert(fabs(B[27] - 0.0807167333268751) < 1e-14,
			"Incorrect value of B at 27");
	mu_assert(fabs(B[28] - 0.0499219869343351) < 1e-14,
			"Incorrect value of B at 28");
	mu_assert(fabs(B[29] - -0.0023736192508975) < 1e-14,
			"Incorrect value of B at 29");

	// end test code //

	free(A);
	free(B);

	return NULL;
}

char *test_dsysv()
{
	int n = 6,
	    m = 5;

	// start test code //
	double *A = Calloc(double, n*n);
	double *B = Calloc(double, n*m);

	// only store the upper triangular part of A
	matrix_set(A, n, 0, 0, 0.4543421836368821);
	matrix_set(A, n, 0, 1, 0.8708338836669620);
	matrix_set(A, n, 0, 2, 1.3638340495356920);
	matrix_set(A, n, 0, 3, 0.8361050302144852);
	matrix_set(A, n, 0, 4, 1.3203463886997013);
	matrix_set(A, n, 0, 5, 0.3915472119381547);
	matrix_set(A, n, 1, 1, 1.4781728513484600);
	matrix_set(A, n, 1, 2, 1.7275151336935415);
	matrix_set(A, n, 1, 3, 1.1817139356024176);
	matrix_set(A, n, 1, 4, 0.7436086782250922);
	matrix_set(A, n, 1, 5, 0.1101758222549450);
	matrix_set(A, n, 2, 2, 1.9363682709237851);
	matrix_set(A, n, 2, 3, 1.1255164391384127);
	matrix_set(A, n, 2, 4, 1.0935575148560115);
	matrix_set(A, n, 2, 5, 1.4678279983625921);
	matrix_set(A, n, 3, 3, 1.7500757162326757);
	matrix_set(A, n, 3, 4, 1.5490921663229316);
	matrix_set(A, n, 3, 5, 1.0305675837706338);
	matrix_set(A, n, 4, 4, 0.4015851628106807);
	matrix_set(A, n, 4, 5, 1.2487496402900566);
	matrix_set(A, n, 5, 5, 0.7245473723012897);

	// Store B in column-major order!
	B[0] = 0.6037912122210694;
	B[1] = 0.5464186020522516;
	B[2] = 0.1810847918541411;
	B[3] = 0.1418268895582175;
	B[4] = 0.5459836535934901;
	B[5] = 0.5890609930309275;
	B[6] = 0.1128454279279324;
	B[7] = 0.8930541056550655;
	B[8] = 0.6946437745982983;
	B[9] = 0.0955380494302154;
	B[10] = 0.5750037200376288;
	B[11] = 0.0326245221201559;
	B[12] = 0.3336701777312929;
	B[13] = 0.7648765739095678;
	B[14] = 0.2662986413718805;
	B[15] = 0.7850810368985298;
	B[16] = 0.5432388739552745;
	B[17] = 0.4387739258059151;
	B[18] = 0.4257906469646436;
	B[19] = 0.1272470768775465;
	B[20] = 0.4276616397814972;
	B[21] = 0.8137579718316245;
	B[22] = 0.6849003723960281;
	B[23] = 0.1768571691078990;
	B[24] = 0.4183278358395650;
	B[25] = 0.6517633972400351;
	B[26] = 0.1154775550239331;
	B[27] = 0.4217248849174023;
	B[28] = 0.9179697263236190;
	B[29] = 0.6532254399609347;

	// run dsysv, note that Lapack expects matrices to be in column-major
	// order, so we can use the 'L' notation for the triangle of A, since
	// we've stored the upper triangle in Row-Major order.

	int *IPIV = Calloc(int, n);
	double *WORK = Calloc(double, 1);
	int status;

	// first we determine the necessary size of the WORK array
	status = dsysv('L', n, m, A, n, IPIV, B, n, WORK, -1);
	mu_assert(status == 0, "dsysv workspace query failed");

	int LWORK = WORK[0];
	WORK = Realloc(WORK, double, LWORK);
	status = dsysv('L', n, m, A, n, IPIV, B, n, WORK, LWORK);
	mu_assert(status == 0, "dsysv didn't return status success");

	// Since B now contains the solution in Column-Major order, we have to
	// check it in that order

	mu_assert(fabs(B[0] - 3.0915448286548806) < 1e-14,
			"Incorrect value of B at 0");
	mu_assert(fabs(B[1] - -1.2114333666218096) < 1e-14,
			"Incorrect value of B at 1");
	mu_assert(fabs(B[2] - -0.1734297056577389) < 1e-14,
			"Incorrect value of B at 2");
	mu_assert(fabs(B[3] - -0.6989941801726605) < 1e-14,
			"Incorrect value of B at 3");
	mu_assert(fabs(B[4] - 1.2577948324106381) < 1e-14,
			"Incorrect value of B at 4");
	mu_assert(fabs(B[5] - -1.4956913279293909) < 1e-14,
			"Incorrect value of B at 5");
	mu_assert(fabs(B[6] - -0.2923881304345451) < 1e-14,
			"Incorrect value of B at 6");
	mu_assert(fabs(B[7] - 0.3467010144627596) < 1e-14,
			"Incorrect value of B at 7");
	mu_assert(fabs(B[8] - 0.4892730831569431) < 1e-14,
			"Incorrect value of B at 8");
	mu_assert(fabs(B[9] - 0.2811039364176572) < 1e-14,
			"Incorrect value of B at 9");
	mu_assert(fabs(B[10] - -0.7323586733947237) < 1e-14,
			"Incorrect value of B at 10");
	mu_assert(fabs(B[11] - 0.0214996365534143) < 1e-14,
			"Incorrect value of B at 11");
	mu_assert(fabs(B[12] - -0.9355264353773129) < 1e-14,
			"Incorrect value of B at 12");
	mu_assert(fabs(B[13] - 0.2709743256273919) < 1e-14,
			"Incorrect value of B at 13");
	mu_assert(fabs(B[14] - 0.2328234557913496) < 1e-14,
			"Incorrect value of B at 14");
	mu_assert(fabs(B[15] - 0.9367092697976086) < 1e-14,
			"Incorrect value of B at 15");
	mu_assert(fabs(B[16] - -0.4501075692310449) < 1e-14,
			"Incorrect value of B at 16");
	mu_assert(fabs(B[17] - 0.0416902255163366) < 1e-14,
			"Incorrect value of B at 17");
	mu_assert(fabs(B[18] - 2.2216982312706905) < 1e-14,
			"Incorrect value of B at 18");
	mu_assert(fabs(B[19] - -0.4820931673893176) < 1e-14,
			"Incorrect value of B at 19");
	mu_assert(fabs(B[20] - -0.8456462251088332) < 1e-14,
			"Incorrect value of B at 20");
	mu_assert(fabs(B[21] - -0.3761761825939751) < 1e-14,
			"Incorrect value of B at 21");
	mu_assert(fabs(B[22] - 1.1921920764994696) < 1e-14,
			"Incorrect value of B at 22");
	mu_assert(fabs(B[23] - -0.6897255145640184) < 1e-14,
			"Incorrect value of B at 23");
	mu_assert(fabs(B[24] - 2.0325624816957180) < 1e-14,
			"Incorrect value of B at 24");
	mu_assert(fabs(B[25] - -0.9900930297944344) < 1e-14,
			"Incorrect value of B at 25");
	mu_assert(fabs(B[26] - -0.0462533459389938) < 1e-14,
			"Incorrect value of B at 26");
	mu_assert(fabs(B[27] - 0.0960916931433909) < 1e-14,
			"Incorrect value of B at 27");
	mu_assert(fabs(B[28] - 0.5805302287627144) < 1e-14,
			"Incorrect value of B at 28");
	mu_assert(fabs(B[29] - -1.0897953557455400) < 1e-14,
			"Incorrect value of B at 29");

	free(WORK);
	free(IPIV);

	// end test code //

	free(A);
	free(B);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();

	mu_run_test(test_gensvm_calculate_omega);
	mu_run_test(test_gensvm_majorize_is_simple);
	mu_run_test(test_gensvm_calculate_ab_non_simple);
	mu_run_test(test_gensvm_calculate_ab_simple);

	mu_run_test(test_dposv);
	mu_run_test(test_dsysv);

	mu_run_test(test_gensvm_get_update);
	mu_run_test(test_gensvm_get_update_sparse);

	return NULL;
}

RUN_TESTS(all_tests);
