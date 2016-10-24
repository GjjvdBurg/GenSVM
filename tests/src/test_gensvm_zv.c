/**
 * @file test_gensvm_zv.c
 * @author G.J.J. van den Burg
 * @date 2016-10-17
 * @brief Unit tests for gensvm_zv.c functions
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
#include "gensvm_zv.h"

char *test_zv_dense_1()
{
	int n = 8,
	    m = 3,
	    K = 3;

	struct GenModel *model = gensvm_init_model();
	model->n = n;
	model->m = m;
	model->K = K;
	model->V = Calloc(double, (m+1)*(K-1));
	matrix_set(model->V, model->K-1, 0, 0, 0.9025324416711976);
	matrix_set(model->V, model->K-1, 0, 1, 0.9776784486541952);
	matrix_set(model->V, model->K-1, 1, 0, 0.8336347240271171);
	matrix_set(model->V, model->K-1, 1, 1, 0.1213543508830703);
	matrix_set(model->V, model->K-1, 2, 0, 0.9401310852208050);
	matrix_set(model->V, model->K-1, 2, 1, 0.7407478086613410);
	matrix_set(model->V, model->K-1, 3, 0, 0.9053353815353901);
	matrix_set(model->V, model->K-1, 3, 1, 0.8056059951641629);

	struct GenData *data = gensvm_init_data();
	data->n = n;
	data->m = m;
	data->K = K;
	data->Z = Calloc(double, n*(m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.4787662921736276);
	matrix_set(data->Z, data->m+1, 0, 2, 0.7983044792882817);
	matrix_set(data->Z, data->m+1, 0, 3, 0.4273006962165122);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.7160319769123790);
	matrix_set(data->Z, data->m+1, 1, 2, 0.5233066338418962);
	matrix_set(data->Z, data->m+1, 1, 3, 0.4063256860579537);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, 0.3735389652435536);
	matrix_set(data->Z, data->m+1, 2, 2, 0.8156214578257802);
	matrix_set(data->Z, data->m+1, 2, 3, 0.6928367712901857);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, 0.3694690105850765);
	matrix_set(data->Z, data->m+1, 3, 2, 0.8539671806454873);
	matrix_set(data->Z, data->m+1, 3, 3, 0.5455108033084728);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, 0.8802158533820680);
	matrix_set(data->Z, data->m+1, 4, 2, 0.0690778177684403);
	matrix_set(data->Z, data->m+1, 4, 3, 0.4513353324958240);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, 0.7752402729955837);
	matrix_set(data->Z, data->m+1, 5, 2, 0.3941285577056867);
	matrix_set(data->Z, data->m+1, 5, 3, 0.2921042477960945);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, 0.6139038657913901);
	matrix_set(data->Z, data->m+1, 6, 2, 0.4529743309354828);
	matrix_set(data->Z, data->m+1, 6, 3, 0.7295983135133345);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, 0.7663625136928905);
	matrix_set(data->Z, data->m+1, 7, 2, 0.3845759571625976);
	matrix_set(data->Z, data->m+1, 7, 3, 0.2291505633226144);

	// start test code //
	double *ZV = Calloc(double, n*(K-1));
	double eps = 1e-14;
	gensvm_calculate_ZV(model, data, ZV);

	mu_assert(fabs(matrix_get(ZV, K-1, 0, 0) - 2.4390099428102818) < eps,
			"Incorrect ZV at 0, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 0, 1) - 1.9713571175527906) < eps,
			"Incorrect ZV at 0, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 0) - 2.3592794147310747) < eps,
			"Incorrect ZV at 1, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 1) - 1.7795486953777246) < eps,
			"Incorrect ZV at 1, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 0) - 2.6079682228282564) < eps,
			"Incorrect ZV at 2, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 1) - 2.1853322915140310) < eps,
			"Incorrect ZV at 2, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 3, 0) - 2.5072459618750060) < eps,
			"Incorrect ZV at 3, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 3, 1) - 2.0945562119091297) < eps,
			"Incorrect ZV at 3, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 4, 0) - 2.1098629909184887) < eps,
			"Incorrect ZV at 4, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 4, 1) - 1.4992641640054902) < eps,
			"Incorrect ZV at 4, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 5, 0) - 2.1837844720035213) < eps,
			"Incorrect ZV at 5, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 5, 1) - 1.5990280274507829) < eps,
			"Incorrect ZV at 5, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 6, 0) - 2.5006904382610986) < eps,
			"Incorrect ZV at 6, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 6, 1) - 1.9754868722402175) < eps,
			"Incorrect ZV at 6, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 7, 0) - 2.1104087689101294) < eps,
			"Incorrect ZV at 7, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 7, 1) - 1.5401587391844891) < eps,
			"Incorrect ZV at 7, 1");

	free(ZV);
	// end test code //
	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_zv_dense_2()
{
	int n = 8,
	    m = 3,
	    K = 3;

	struct GenModel *model = gensvm_init_model();
	model->n = n;
	model->m = m;
	model->K = K;
	model->V = Calloc(double, (m+1)*(K-1));
	matrix_set(model->V, model->K-1, 0, 0, 0.9025324416711976);
	matrix_set(model->V, model->K-1, 0, 1, 0.9776784486541952);
	matrix_set(model->V, model->K-1, 1, 0, 0.8336347240271171);
	matrix_set(model->V, model->K-1, 1, 1, 0.1213543508830703);
	matrix_set(model->V, model->K-1, 2, 0, 0.9401310852208050);
	matrix_set(model->V, model->K-1, 2, 1, 0.7407478086613410);
	matrix_set(model->V, model->K-1, 3, 0, 0.9053353815353901);
	matrix_set(model->V, model->K-1, 3, 1, 0.8056059951641629);

	struct GenData *data = gensvm_init_data();
	data->n = 3;
	data->m = m;
	data->K = K;
	data->Z = Calloc(double, data->n*(data->m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.4787662921736276);
	matrix_set(data->Z, data->m+1, 0, 2, 0.7983044792882817);
	matrix_set(data->Z, data->m+1, 0, 3, 0.4273006962165122);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.7160319769123790);
	matrix_set(data->Z, data->m+1, 1, 2, 0.5233066338418962);
	matrix_set(data->Z, data->m+1, 1, 3, 0.4063256860579537);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, 0.3735389652435536);
	matrix_set(data->Z, data->m+1, 2, 2, 0.8156214578257802);
	matrix_set(data->Z, data->m+1, 2, 3, 0.6928367712901857);

	// start test code //
	double *ZV = Calloc(double, data->n*(K-1));
	double eps = 1e-14;
	gensvm_calculate_ZV(model, data, ZV);

	mu_assert(fabs(matrix_get(ZV, K-1, 0, 0) - 2.4390099428102818) < eps,
			"Incorrect ZV at 0, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 0, 1) - 1.9713571175527906) < eps,
			"Incorrect ZV at 0, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 0) - 2.3592794147310747) < eps,
			"Incorrect ZV at 1, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 1) - 1.7795486953777246) < eps,
			"Incorrect ZV at 1, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 0) - 2.6079682228282564) < eps,
			"Incorrect ZV at 2, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 1) - 2.1853322915140310) < eps,
			"Incorrect ZV at 2, 1");

	free(ZV);
	// end test code //
	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_zv_sparse_1()
{
	int n = 8,
	    m = 3,
	    K = 3;

	struct GenModel *model = gensvm_init_model();
	model->n = n;
	model->m = m;
	model->K = K;
	model->V = Calloc(double, (m+1)*(K-1));
	matrix_set(model->V, model->K-1, 0, 0, 0.9025324416711976);
	matrix_set(model->V, model->K-1, 0, 1, 0.9776784486541952);
	matrix_set(model->V, model->K-1, 1, 0, 0.8336347240271171);
	matrix_set(model->V, model->K-1, 1, 1, 0.1213543508830703);
	matrix_set(model->V, model->K-1, 2, 0, 0.9401310852208050);
	matrix_set(model->V, model->K-1, 2, 1, 0.7407478086613410);
	matrix_set(model->V, model->K-1, 3, 0, 0.9053353815353901);
	matrix_set(model->V, model->K-1, 3, 1, 0.8056059951641629);

	struct GenData *data = gensvm_init_data();
	data->n = n;
	data->m = m;
	data->K = K;
	data->Z = Calloc(double, n*(m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.4787662921736276);
	matrix_set(data->Z, data->m+1, 0, 2, 0.7983044792882817);
	matrix_set(data->Z, data->m+1, 0, 3, 0.4273006962165122);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.7160319769123790);
	matrix_set(data->Z, data->m+1, 1, 2, 0.5233066338418962);
	matrix_set(data->Z, data->m+1, 1, 3, 0.4063256860579537);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, 0.3735389652435536);
	matrix_set(data->Z, data->m+1, 2, 2, 0.8156214578257802);
	matrix_set(data->Z, data->m+1, 2, 3, 0.6928367712901857);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, 0.3694690105850765);
	matrix_set(data->Z, data->m+1, 3, 2, 0.8539671806454873);
	matrix_set(data->Z, data->m+1, 3, 3, 0.5455108033084728);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, 0.8802158533820680);
	matrix_set(data->Z, data->m+1, 4, 2, 0.0690778177684403);
	matrix_set(data->Z, data->m+1, 4, 3, 0.4513353324958240);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, 0.7752402729955837);
	matrix_set(data->Z, data->m+1, 5, 2, 0.3941285577056867);
	matrix_set(data->Z, data->m+1, 5, 3, 0.2921042477960945);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, 0.6139038657913901);
	matrix_set(data->Z, data->m+1, 6, 2, 0.4529743309354828);
	matrix_set(data->Z, data->m+1, 6, 3, 0.7295983135133345);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, 0.7663625136928905);
	matrix_set(data->Z, data->m+1, 7, 2, 0.3845759571625976);
	matrix_set(data->Z, data->m+1, 7, 3, 0.2291505633226144);

	// convert Z to sparse matrix
	data->spZ = gensvm_dense_to_sparse(data->Z, data->n, data->m+1);
	free(data->Z);
	data->RAW = NULL;
	data->Z = NULL;

	// start test code //
	double *ZV = Calloc(double, n*(K-1));
	double eps = 1e-14;
	gensvm_calculate_ZV(model, data, ZV);

	mu_assert(fabs(matrix_get(ZV, K-1, 0, 0) - 2.4390099428102818) < eps,
			"Incorrect ZV at 0, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 0, 1) - 1.9713571175527906) < eps,
			"Incorrect ZV at 0, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 0) - 2.3592794147310747) < eps,
			"Incorrect ZV at 1, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 1) - 1.7795486953777246) < eps,
			"Incorrect ZV at 1, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 0) - 2.6079682228282564) < eps,
			"Incorrect ZV at 2, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 1) - 2.1853322915140310) < eps,
			"Incorrect ZV at 2, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 3, 0) - 2.5072459618750060) < eps,
			"Incorrect ZV at 3, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 3, 1) - 2.0945562119091297) < eps,
			"Incorrect ZV at 3, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 4, 0) - 2.1098629909184887) < eps,
			"Incorrect ZV at 4, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 4, 1) - 1.4992641640054902) < eps,
			"Incorrect ZV at 4, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 5, 0) - 2.1837844720035213) < eps,
			"Incorrect ZV at 5, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 5, 1) - 1.5990280274507829) < eps,
			"Incorrect ZV at 5, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 6, 0) - 2.5006904382610986) < eps,
			"Incorrect ZV at 6, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 6, 1) - 1.9754868722402175) < eps,
			"Incorrect ZV at 6, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 7, 0) - 2.1104087689101294) < eps,
			"Incorrect ZV at 7, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 7, 1) - 1.5401587391844891) < eps,
			"Incorrect ZV at 7, 1");

	free(ZV);
	// end test code //
	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_zv_sparse_2()
{
	int n = 8,
	    m = 3,
	    K = 3;

	struct GenModel *model = gensvm_init_model();
	model->n = n;
	model->m = m;
	model->K = K;
	model->V = Calloc(double, (m+1)*(K-1));
	matrix_set(model->V, model->K-1, 0, 0, 0.9025324416711976);
	matrix_set(model->V, model->K-1, 0, 1, 0.9776784486541952);
	matrix_set(model->V, model->K-1, 1, 0, 0.8336347240271171);
	matrix_set(model->V, model->K-1, 1, 1, 0.1213543508830703);
	matrix_set(model->V, model->K-1, 2, 0, 0.9401310852208050);
	matrix_set(model->V, model->K-1, 2, 1, 0.7407478086613410);
	matrix_set(model->V, model->K-1, 3, 0, 0.9053353815353901);
	matrix_set(model->V, model->K-1, 3, 1, 0.8056059951641629);

	struct GenData *data = gensvm_init_data();
	data->n = 3;
	data->m = m;
	data->K = K;
	data->Z = Calloc(double, data->n*(data->m+1));
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.4787662921736276);
	matrix_set(data->Z, data->m+1, 0, 2, 0.7983044792882817);
	matrix_set(data->Z, data->m+1, 0, 3, 0.4273006962165122);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.7160319769123790);
	matrix_set(data->Z, data->m+1, 1, 2, 0.5233066338418962);
	matrix_set(data->Z, data->m+1, 1, 3, 0.4063256860579537);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, 0.3735389652435536);
	matrix_set(data->Z, data->m+1, 2, 2, 0.8156214578257802);
	matrix_set(data->Z, data->m+1, 2, 3, 0.6928367712901857);

	// convert Z to sparse matrix
	data->spZ = gensvm_dense_to_sparse(data->Z, data->n, data->m+1);
	free(data->Z);
	data->RAW = NULL;
	data->Z = NULL;

	// start test code //
	double *ZV = Calloc(double, data->n*(data->K-1));
	double eps = 1e-14;
	gensvm_calculate_ZV(model, data, ZV);

	mu_assert(fabs(matrix_get(ZV, K-1, 0, 0) - 2.4390099428102818) < eps,
			"Incorrect ZV at 0, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 0, 1) - 1.9713571175527906) < eps,
			"Incorrect ZV at 0, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 0) - 2.3592794147310747) < eps,
			"Incorrect ZV at 1, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 1, 1) - 1.7795486953777246) < eps,
			"Incorrect ZV at 1, 1");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 0) - 2.6079682228282564) < eps,
			"Incorrect ZV at 2, 0");
	mu_assert(fabs(matrix_get(ZV, K-1, 2, 1) - 2.1853322915140310) < eps,
			"Incorrect ZV at 2, 1");

	free(ZV);
	// end test code //
	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_zv_dense_1);
	mu_run_test(test_zv_dense_2);
	mu_run_test(test_zv_sparse_1);
	mu_run_test(test_zv_sparse_2);

	return NULL;
}

RUN_TESTS(all_tests);
