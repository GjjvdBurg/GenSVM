/**
 * @file test_gensvm_optimize.c
 * @author Gertjan van den Burg
 * @date September, 2016
 * @brief Unit tests for gensvm_optimize.c functions
 */

#include "minunit.h"
#include "gensvm_optimize.h"

char *test_gensvm_optimize()
{
	mu_test_missing();
	return NULL;
}

char *test_gensvm_get_loss()
{
	mu_test_missing();
	return NULL;
}

char *test_gensvm_get_update()
{
	mu_test_missing();
	return NULL;
}

char *test_gensvm_category_matrix()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	int n = 5,
	    m = 3,
	    K = 3;

	data->n = n;
	data->m = m;
	data->K = K;

	model->n = n;
	model->m = m;
	model->K = K;

	gensvm_allocate_model(model);
	data->y = Calloc(long, data->n);
	data->y[0] = 1;
	data->y[1] = 2;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 1;

	// start test code //

	gensvm_category_matrix(model, data);

	mu_assert(matrix_get(model->R, K, 0, 0) == 0, "Incorrect R at 0, 0");
	mu_assert(matrix_get(model->R, K, 0, 1) == 1, "Incorrect R at 0, 1");
	mu_assert(matrix_get(model->R, K, 0, 2) == 1, "Incorrect R at 0, 2");

	mu_assert(matrix_get(model->R, K, 1, 0) == 1, "Incorrect R at 1, 0");
	mu_assert(matrix_get(model->R, K, 1, 1) == 0, "Incorrect R at 1, 1");
	mu_assert(matrix_get(model->R, K, 1, 2) == 1, "Incorrect R at 1, 2");

	mu_assert(matrix_get(model->R, K, 2, 0) == 1, "Incorrect R at 2, 0");
	mu_assert(matrix_get(model->R, K, 2, 1) == 1, "Incorrect R at 2, 1");
	mu_assert(matrix_get(model->R, K, 2, 2) == 0, "Incorrect R at 2, 2");

	mu_assert(matrix_get(model->R, K, 3, 0) == 1, "Incorrect R at 3, 0");
	mu_assert(matrix_get(model->R, K, 3, 1) == 0, "Incorrect R at 3, 1");
	mu_assert(matrix_get(model->R, K, 3, 2) == 1, "Incorrect R at 3, 2");

	mu_assert(matrix_get(model->R, K, 4, 0) == 0, "Incorrect R at 4, 0");
	mu_assert(matrix_get(model->R, K, 4, 1) == 1, "Incorrect R at 4, 1");
	mu_assert(matrix_get(model->R, K, 4, 2) == 1, "Incorrect R at 4, 2");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_simplex_diff()
{
	struct GenData *data = gensvm_init_data();
	struct GenModel *model = gensvm_init_model();

	int n = 8,
	    m = 3,
	    K = 3;
	model->n = n;
	model->m = m;
	model->K = K;
	data->n = n;
	data->m = m;
	data->K = K;

	data->y = Calloc(long, n);

	gensvm_allocate_model(model);
	gensvm_simplex(model->K, model->U);

	data->y[0] = 2;
	data->y[1] = 1;
	data->y[2] = 3;
	data->y[3] = 2;
	data->y[4] = 3;
	data->y[5] = 3;
	data->y[6] = 1;
	data->y[7] = 2;

	// start test code //
	gensvm_simplex_diff(model, data);

	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 0, 0, 0) -
				1.0000000000000000) < 1e-14,
			"Incorrect value at UU(0, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 0, 0, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(0, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 0, 0, 2) -
				0.5000000000000000) < 1e-14,
			"Incorrect value at UU(0, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 0, 1, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(0, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 0, 1, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(0, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 0, 1, 2) -
				-0.8660254037844388) < 1e-14,
			"Incorrect value at UU(0, 1, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 1, 0, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(1, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 1, 0, 1) -
				-1.0000000000000000) < 1e-14,
			"Incorrect value at UU(1, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 1, 0, 2) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value at UU(1, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 1, 1, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(1, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 1, 1, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(1, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 1, 1, 2) -
				-0.8660254037844388) < 1e-14,
			"Incorrect value at UU(1, 1, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 2, 0, 0) -
				0.5000000000000000) < 1e-14,
			"Incorrect value at UU(2, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 2, 0, 1) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value at UU(2, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 2, 0, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(2, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 2, 1, 0) -
				0.8660254037844388) < 1e-14,
			"Incorrect value at UU(2, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 2, 1, 1) -
				0.8660254037844388) < 1e-14,
			"Incorrect value at UU(2, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 2, 1, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(2, 1, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 3, 0, 0) -
				1.0000000000000000) < 1e-14,
			"Incorrect value at UU(3, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 3, 0, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(3, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 3, 0, 2) -
				0.5000000000000000) < 1e-14,
			"Incorrect value at UU(3, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 3, 1, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(3, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 3, 1, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(3, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 3, 1, 2) -
				-0.8660254037844388) < 1e-14,
			"Incorrect value at UU(3, 1, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 4, 0, 0) -
				0.5000000000000000) < 1e-14,
			"Incorrect value at UU(4, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 4, 0, 1) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value at UU(4, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 4, 0, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(4, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 4, 1, 0) -
				0.8660254037844388) < 1e-14,
			"Incorrect value at UU(4, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 4, 1, 1) -
				0.8660254037844388) < 1e-14,
			"Incorrect value at UU(4, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 4, 1, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(4, 1, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 5, 0, 0) -
				0.5000000000000000) < 1e-14,
			"Incorrect value at UU(5, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 5, 0, 1) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value at UU(5, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 5, 0, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(5, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 5, 1, 0) -
				0.8660254037844388) < 1e-14,
			"Incorrect value at UU(5, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 5, 1, 1) -
				0.8660254037844388) < 1e-14,
			"Incorrect value at UU(5, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 5, 1, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(5, 1, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 6, 0, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(6, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 6, 0, 1) -
				-1.0000000000000000) < 1e-14,
			"Incorrect value at UU(6, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 6, 0, 2) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value at UU(6, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 6, 1, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(6, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 6, 1, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(6, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 6, 1, 2) -
				-0.8660254037844388) < 1e-14,
			"Incorrect value at UU(6, 1, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 7, 0, 0) -
				1.0000000000000000) < 1e-14,
			"Incorrect value at UU(7, 0, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 7, 0, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(7, 0, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 7, 0, 2) -
				0.5000000000000000) < 1e-14,
			"Incorrect value at UU(7, 0, 2)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 7, 1, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(7, 1, 0)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 7, 1, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value at UU(7, 1, 1)");
	mu_assert(fabs(matrix3_get(model->UU, K-1, K, 7, 1, 2) -
				-0.8660254037844388) < 1e-14,
			"Incorrect value at UU(7, 1, 2)");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_calculate_errors()
{
	mu_test_missing();
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

char *test_dposv()
{
	mu_test_missing();
	return NULL;
}

char *test_dsysv()
{
	mu_test_missing();
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_gensvm_optimize);
	mu_run_test(test_gensvm_get_loss);
	mu_run_test(test_gensvm_get_update);
	mu_run_test(test_gensvm_category_matrix);
	mu_run_test(test_gensvm_simplex_diff);
	mu_run_test(test_gensvm_calculate_errors);
	mu_run_test(test_gensvm_calculate_huber);
	mu_run_test(test_gensvm_step_doubling);
	mu_run_test(test_dposv);
	mu_run_test(test_dsysv);

	return NULL;
}

RUN_TESTS(all_tests);
