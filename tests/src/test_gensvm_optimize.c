/**
 * @file test_gensvm_optimize.c
 * @author Gertjan van den Burg
 * @date September, 2016
 * @brief Unit tests for gensvm_optimize.c functions
 */

#include "minunit.h"
#include "gensvm_optimize.h"
#include "gensvm_init.h"

char *test_gensvm_optimize()
{
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
	model->weight_idx = 1;
	model->kappa = 0.5;
	model->p = 1.5;
	model->lambda = 0.123;

	gensvm_allocate_model(model);
	gensvm_initialize_weights(data, model);
	gensvm_simplex(model->K, model->U);
	gensvm_simplex_diff(model, data);
	gensvm_category_matrix(model, data);

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
	double loss = gensvm_get_loss(model, data, ZV);

	mu_assert(fabs(loss - 0.903071383013108) < 1e-14,
			"Incorrect value of the loss");

	free(ZV);
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);


	return NULL;
}

char *test_gensvm_get_loss_2()
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
	model->weight_idx = 2;
	model->kappa = 0.5;
	model->p = 1.5;
	model->lambda = 0.123;

	gensvm_allocate_model(model);
	gensvm_initialize_weights(data, model);
	gensvm_simplex(model->K, model->U);
	gensvm_simplex_diff(model, data);
	gensvm_category_matrix(model, data);

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
	double loss = gensvm_get_loss(model, data, ZV);

	mu_assert(fabs(loss - 0.972847045993281) < 1e-14,
			"Incorrect value of the loss");

	free(ZV);
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);


	return NULL;
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
	gensvm_simplex(model->K, model->U);
	gensvm_simplex_diff(model, data);

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
	mu_run_test(test_gensvm_optimize);
	mu_run_test(test_gensvm_get_loss_1);
	mu_run_test(test_gensvm_get_loss_2);
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
