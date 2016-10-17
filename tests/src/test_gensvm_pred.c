/**
 * @file test_gensvm_pred.c
 * @author Gertjan van den Burg
 * @date September, 2016
 * @brief Unit tests for gensvm_pred.c functions
 */
#include "minunit.h"
#include "gensvm_pred.h"

/**
 * This testcase is designed as follows: 12 evenly spaced points are plotted 
 * on the unit circle in the simplex space. These points are the ones for 
 * which we want to predict the class. To get these points, we need a Z and a 
 * V which map to these points. To get this, Z was equal to [1 Q] and V was 
 * equal to [0; R] where Q and R are from the reduced QR decomposition of the 
 * 12x2 matrix S which contains the points in simplex space. Here's the 
 * Matlab/Octave code to generate this data:
 *
 *     n = 12;
 *     K = 3;
 *     S = [cos(1/12*pi+1/6*pi*[0:(n-1)])', sin(1/12*pi+1/6*pi*[0:(n-1)])'];
 *     [Q, R] = qr(S, '0');
 *     Z = [ones(n, 1), Q];
 *     V = [zeros(1, K-1); R];
 *
 */
char *test_gensvm_predict_labels_dense()
{
	int n = 12;
	int m = 2;
	int K = 3;

	struct GenData *data = gensvm_init_data();
	struct GenModel *model = gensvm_init_model();

	model->n = n;
	model->m = m;
	model->K = K;

	data->n = n;
	data->m = m;
	data->r = m;
	data->K = K;

	data->Z = Calloc(double, n*(m+1));
	data->y = Calloc(long, n);

	matrix_set(data->Z, m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 0, 1, -0.3943375672974065);
	matrix_set(data->Z, m+1, 0, 2, -0.1056624327025935);
	matrix_set(data->Z, m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 1, 1, -0.2886751345948129);
	matrix_set(data->Z, m+1, 1, 2, -0.2886751345948128);
	matrix_set(data->Z, m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 2, 1, -0.1056624327025937);
	matrix_set(data->Z, m+1, 2, 2, -0.3943375672974063);
	matrix_set(data->Z, m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 3, 1, 0.1056624327025935);
	matrix_set(data->Z, m+1, 3, 2, -0.3943375672974064);
	matrix_set(data->Z, m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 4, 1, 0.2886751345948129);
	matrix_set(data->Z, m+1, 4, 2, -0.2886751345948129);
	matrix_set(data->Z, m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 5, 1, 0.3943375672974064);
	matrix_set(data->Z, m+1, 5, 2, -0.1056624327025937);
	matrix_set(data->Z, m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 6, 1, 0.3943375672974065);
	matrix_set(data->Z, m+1, 6, 2, 0.1056624327025935);
	matrix_set(data->Z, m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 7, 1, 0.2886751345948130);
	matrix_set(data->Z, m+1, 7, 2, 0.2886751345948128);
	matrix_set(data->Z, m+1, 8, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 8, 1, 0.1056624327025939);
	matrix_set(data->Z, m+1, 8, 2, 0.3943375672974063);
	matrix_set(data->Z, m+1, 9, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 9, 1, -0.1056624327025934);
	matrix_set(data->Z, m+1, 9, 2, 0.3943375672974064);
	matrix_set(data->Z, m+1, 10, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 10, 1, -0.2886751345948126);
	matrix_set(data->Z, m+1, 10, 2, 0.2886751345948132);
	matrix_set(data->Z, m+1, 11, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 11, 1, -0.3943375672974064);
	matrix_set(data->Z, m+1, 11, 2, 0.1056624327025939);

	gensvm_allocate_model(model);

	matrix_set(model->V, K-1, 0, 0, 0.0000000000000000);
	matrix_set(model->V, K-1, 0, 1, 0.0000000000000000);
	matrix_set(model->V, K-1, 1, 0, -2.4494897427831779);
	matrix_set(model->V, K-1, 1, 1, -0.0000000000000002);
	matrix_set(model->V, K-1, 2, 0, 0.0000000000000000);
	matrix_set(model->V, K-1, 2, 1, -2.4494897427831783);

	// start test code
	long *predy = Calloc(long, n);
	gensvm_predict_labels(data, model, predy);

	mu_assert(predy[0] == 2, "Incorrect label at index 0");
	mu_assert(predy[1] == 3, "Incorrect label at index 1");
	mu_assert(predy[2] == 3, "Incorrect label at index 2");
	mu_assert(predy[3] == 3, "Incorrect label at index 3");
	mu_assert(predy[4] == 3, "Incorrect label at index 4");
	mu_assert(predy[5] == 1, "Incorrect label at index 5");
	mu_assert(predy[6] == 1, "Incorrect label at index 6");
	mu_assert(predy[7] == 1, "Incorrect label at index 7");
	mu_assert(predy[8] == 1, "Incorrect label at index 8");
	mu_assert(predy[9] == 2, "Incorrect label at index 9");
	mu_assert(predy[10] == 2, "Incorrect label at index 10");
	mu_assert(predy[11] == 2, "Incorrect label at index 11");

	// end test code
	gensvm_free_data(data);
	gensvm_free_model(model);
	free(predy);

	return NULL;
}

char *test_gensvm_predict_labels_sparse()
{
	int n = 12;
	int m = 2;
	int K = 3;

	struct GenData *data = gensvm_init_data();
	struct GenModel *model = gensvm_init_model();

	model->n = n;
	model->m = m;
	model->K = K;

	data->n = n;
	data->m = m;
	data->r = m;
	data->K = K;

	data->Z = Calloc(double, n*(m+1));
	data->y = Calloc(long, n);

	matrix_set(data->Z, m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 0, 1, -0.3943375672974065);
	matrix_set(data->Z, m+1, 0, 2, -0.1056624327025935);
	matrix_set(data->Z, m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 1, 1, -0.2886751345948129);
	matrix_set(data->Z, m+1, 1, 2, -0.2886751345948128);
	matrix_set(data->Z, m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 2, 1, -0.1056624327025937);
	matrix_set(data->Z, m+1, 2, 2, -0.3943375672974063);
	matrix_set(data->Z, m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 3, 1, 0.1056624327025935);
	matrix_set(data->Z, m+1, 3, 2, -0.3943375672974064);
	matrix_set(data->Z, m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 4, 1, 0.2886751345948129);
	matrix_set(data->Z, m+1, 4, 2, -0.2886751345948129);
	matrix_set(data->Z, m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 5, 1, 0.3943375672974064);
	matrix_set(data->Z, m+1, 5, 2, -0.1056624327025937);
	matrix_set(data->Z, m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 6, 1, 0.3943375672974065);
	matrix_set(data->Z, m+1, 6, 2, 0.1056624327025935);
	matrix_set(data->Z, m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 7, 1, 0.2886751345948130);
	matrix_set(data->Z, m+1, 7, 2, 0.2886751345948128);
	matrix_set(data->Z, m+1, 8, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 8, 1, 0.1056624327025939);
	matrix_set(data->Z, m+1, 8, 2, 0.3943375672974063);
	matrix_set(data->Z, m+1, 9, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 9, 1, -0.1056624327025934);
	matrix_set(data->Z, m+1, 9, 2, 0.3943375672974064);
	matrix_set(data->Z, m+1, 10, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 10, 1, -0.2886751345948126);
	matrix_set(data->Z, m+1, 10, 2, 0.2886751345948132);
	matrix_set(data->Z, m+1, 11, 0, 1.0000000000000000);
	matrix_set(data->Z, m+1, 11, 1, -0.3943375672974064);
	matrix_set(data->Z, m+1, 11, 2, 0.1056624327025939);

	// convert Z to a sparse matrix to test the sparse functions
	data->spZ = gensvm_dense_to_sparse(data->Z, data->n, data->m+1);
	free(data->Z);
	data->RAW = NULL;
	data->Z = NULL;

	gensvm_allocate_model(model);

	matrix_set(model->V, K-1, 0, 0, 0.0000000000000000);
	matrix_set(model->V, K-1, 0, 1, 0.0000000000000000);
	matrix_set(model->V, K-1, 1, 0, -2.4494897427831779);
	matrix_set(model->V, K-1, 1, 1, -0.0000000000000002);
	matrix_set(model->V, K-1, 2, 0, 0.0000000000000000);
	matrix_set(model->V, K-1, 2, 1, -2.4494897427831783);

	// start test code
	long *predy = Calloc(long, n);
	gensvm_predict_labels(data, model, predy);

	mu_assert(predy[0] == 2, "Incorrect label at index 0");
	mu_assert(predy[1] == 3, "Incorrect label at index 1");
	mu_assert(predy[2] == 3, "Incorrect label at index 2");
	mu_assert(predy[3] == 3, "Incorrect label at index 3");
	mu_assert(predy[4] == 3, "Incorrect label at index 4");
	mu_assert(predy[5] == 1, "Incorrect label at index 5");
	mu_assert(predy[6] == 1, "Incorrect label at index 6");
	mu_assert(predy[7] == 1, "Incorrect label at index 7");
	mu_assert(predy[8] == 1, "Incorrect label at index 8");
	mu_assert(predy[9] == 2, "Incorrect label at index 9");
	mu_assert(predy[10] == 2, "Incorrect label at index 10");
	mu_assert(predy[11] == 2, "Incorrect label at index 11");

	// end test code
	gensvm_free_data(data);
	gensvm_free_model(model);
	free(predy);

	return NULL;
}

char *test_gensvm_prediction_perf()
{
	int i, n = 8;
	struct GenData *data = gensvm_init_data();
	data->n = n;
	data->y = Calloc(long, n);
	data->y[0] = 1;
	data->y[1] = 1;
	data->y[2] = 1;
	data->y[3] = 1;
	data->y[4] = 2;
	data->y[5] = 2;
	data->y[6] = 2;
	data->y[7] = 3;

	long *y = Calloc(long, n);
	for (i=0; i<n; i++)
		y[i] = 1;
	mu_assert(gensvm_prediction_perf(data, y) == 50.0,
			"Incorrect first time.");

	for (i=0; i<n; i++)
		y[i] = 2;
	mu_assert(gensvm_prediction_perf(data, y) == 37.5,
			"Incorrect second time.");

	for (i=0; i<n; i++)
		y[i] = 3;
	mu_assert(gensvm_prediction_perf(data, y) == 12.5,
			"Incorrect third time.");

	free(y);
	gensvm_free_data(data);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_gensvm_predict_labels_dense);
	mu_run_test(test_gensvm_predict_labels_sparse);
	mu_run_test(test_gensvm_prediction_perf);

	return NULL;
}

RUN_TESTS(all_tests);
