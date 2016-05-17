/**
 * @file test_gensvm_copy.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_copy.c functions
 */

#include "minunit.h"
#include "gensvm_copy.h"

char *test_copy_model_linear()
{
	struct GenModel *from_model = gensvm_init_model();
	struct GenModel *to_model = gensvm_init_model();

	from_model->p = 2.0;
	from_model->lambda = 1.0;
	from_model->epsilon = 1e-8;
	from_model->kappa = 1.0;
	from_model->weight_idx = 2;
	from_model->kerneltype = K_LINEAR;

	gensvm_copy_model(from_model, to_model);

	mu_assert(to_model->p == 2.0, "to_model->p incorrect.");
	mu_assert(to_model->lambda == 1.0, "to_model->lambda incorrect.");
	mu_assert(to_model->epsilon == 1e-8, "to_model->epsilon incorrect.");
	mu_assert(to_model->kappa == 1.0, "to_model->kappa incorrect.");
	mu_assert(to_model->weight_idx == 2,
		       	"to_model->weight_idx incorrect.");
	mu_assert(to_model->kerneltype == K_LINEAR, "to->kerneltype incorrect");

	gensvm_free_model(from_model);
	gensvm_free_model(to_model);

	return NULL;
}

char *test_copy_model_poly_1()
{
	struct GenModel *from_model = gensvm_init_model();
	struct GenModel *to_model = gensvm_init_model();

	from_model->p = 2.0;
	from_model->lambda = 1.0;
	from_model->epsilon = 1e-8;
	from_model->kappa = 1.0;
	from_model->weight_idx = 2;
	from_model->kerneltype = K_POLY;
	from_model->kernelparam = Calloc(double, 3);
	from_model->kernelparam[0] = 1.0;
	from_model->kernelparam[1] = 2.0;
	from_model->kernelparam[2] = 3.0;

	gensvm_copy_model(from_model, to_model);

	mu_assert(to_model->p == 2.0, "to->p incorrect.");
	mu_assert(to_model->lambda == 1.0, "to->lambda incorrect.");
	mu_assert(to_model->epsilon == 1e-8, "to->epsilon incorrect.");
	mu_assert(to_model->kappa == 1.0, "to->kappa incorrect.");
	mu_assert(to_model->weight_idx == 2, "to->weight_idx incorrect.");
	mu_assert(to_model->kerneltype == K_POLY, "to->kerneltype incorrect");
	mu_assert(to_model->kernelparam[0] == 1.0,
			"to->kernelparam[0] is incorrect.");
	mu_assert(to_model->kernelparam[1] == 2.0,
			"to->kernelparam[1] is incorrect.");
	mu_assert(to_model->kernelparam[2] == 3.0,
			"to->kernelparam[2] is incorrect.");

	gensvm_free_model(from_model);
	gensvm_free_model(to_model);

	return NULL;
}

char *test_copy_model_poly_2()
{
	struct GenModel *from_model = gensvm_init_model();
	struct GenModel *to_model = gensvm_init_model();

	from_model->p = 2.0;
	from_model->lambda = 1.0;
	from_model->epsilon = 1e-8;
	from_model->kappa = 1.0;
	from_model->weight_idx = 2;
	from_model->kerneltype = K_POLY;
	from_model->kernelparam = Calloc(double, 3);
	from_model->kernelparam[0] = 1.0;
	from_model->kernelparam[1] = 2.0;
	from_model->kernelparam[2] = 3.0;
	to_model->kernelparam = Calloc(double, 3);

	gensvm_copy_model(from_model, to_model);

	mu_assert(to_model->p == 2.0, "to_model->p incorrect.");
	mu_assert(to_model->lambda == 1.0, "to_model->lambda incorrect.");
	mu_assert(to_model->epsilon == 1e-8, "to_model->epsilon incorrect.");
	mu_assert(to_model->kappa == 1.0, "to_model->kappa incorrect.");
	mu_assert(to_model->weight_idx == 2,
		       	"to_model->weight_idx incorrect.");
	mu_assert(to_model->kerneltype == K_POLY, "to->kerneltype incorrect");
	mu_assert(to_model->kernelparam[0] == 1.0,
			"to->kernelparam[0] is incorrect.");
	mu_assert(to_model->kernelparam[1] == 2.0,
			"to->kernelparam[1] is incorrect.");
	mu_assert(to_model->kernelparam[2] == 3.0,
			"to->kernelparam[2] is incorrect.");

	gensvm_free_model(from_model);
	gensvm_free_model(to_model);

	return NULL;
}

char *test_copy_model_rbf_1()
{
	struct GenModel *from_model = gensvm_init_model();
	struct GenModel *to_model = gensvm_init_model();

	from_model->p = 2.0;
	from_model->lambda = 1.0;
	from_model->epsilon = 1e-8;
	from_model->kappa = 1.0;
	from_model->weight_idx = 2;
	from_model->kerneltype = K_RBF;
	from_model->kernelparam = Calloc(double, 1);
	from_model->kernelparam[0] = 5.0;

	gensvm_copy_model(from_model, to_model);

	mu_assert(to_model->p == 2.0, "to_model->p incorrect.");
	mu_assert(to_model->lambda == 1.0, "to_model->lambda incorrect.");
	mu_assert(to_model->epsilon == 1e-8, "to_model->epsilon incorrect.");
	mu_assert(to_model->kappa == 1.0, "to_model->kappa incorrect.");
	mu_assert(to_model->weight_idx == 2,
		       	"to_model->weight_idx incorrect.");
	mu_assert(to_model->kerneltype == K_RBF, "to->kerneltype incorrect");
	mu_assert(to_model->kernelparam[0] == 5.0,
			"to->kernelparam[0] is incorrect.");

	gensvm_free_model(from_model);
	gensvm_free_model(to_model);

	return NULL;
}

char *test_copy_model_rbf_2()
{
	struct GenModel *from_model = gensvm_init_model();
	struct GenModel *to_model = gensvm_init_model();

	from_model->p = 2.0;
	from_model->lambda = 1.0;
	from_model->epsilon = 1e-8;
	from_model->kappa = 1.0;
	from_model->weight_idx = 2;
	from_model->kerneltype = K_RBF;
	from_model->kernelparam = Calloc(double, 1);
	from_model->kernelparam[0] = 5.0;
	to_model->kernelparam = Calloc(double, 1);

	gensvm_copy_model(from_model, to_model);

	mu_assert(to_model->p == 2.0, "to_model->p incorrect.");
	mu_assert(to_model->lambda == 1.0, "to_model->lambda incorrect.");
	mu_assert(to_model->epsilon == 1e-8, "to_model->epsilon incorrect.");
	mu_assert(to_model->kappa == 1.0, "to_model->kappa incorrect.");
	mu_assert(to_model->weight_idx == 2,
		       	"to_model->weight_idx incorrect.");
	mu_assert(to_model->kerneltype == K_RBF, "to->kerneltype incorrect");
	mu_assert(to_model->kernelparam[0] == 5.0,
			"to->kernelparam[0] is incorrect.");

	gensvm_free_model(from_model);
	gensvm_free_model(to_model);

	return NULL;
}

char *test_copy_model_sigmoid_1()
{
	struct GenModel *from_model = gensvm_init_model();
	struct GenModel *to_model = gensvm_init_model();

	from_model->p = 2.0;
	from_model->lambda = 1.0;
	from_model->epsilon = 1e-8;
	from_model->kappa = 1.0;
	from_model->weight_idx = 2;
	from_model->kerneltype = K_SIGMOID;
	from_model->kernelparam = Calloc(double, 2);
	from_model->kernelparam[0] = 5.0;
	from_model->kernelparam[1] = 10.0;

	gensvm_copy_model(from_model, to_model);

	mu_assert(to_model->p == 2.0, "to_model->p incorrect.");
	mu_assert(to_model->lambda == 1.0, "to_model->lambda incorrect.");
	mu_assert(to_model->epsilon == 1e-8, "to_model->epsilon incorrect.");
	mu_assert(to_model->kappa == 1.0, "to_model->kappa incorrect.");
	mu_assert(to_model->weight_idx == 2,
		       	"to_model->weight_idx incorrect.");
	mu_assert(to_model->kerneltype == K_SIGMOID,
		       	"to->kerneltype incorrect");
	mu_assert(to_model->kernelparam[0] == 5.0,
			"to->kernelparam[0] is incorrect.");
	mu_assert(to_model->kernelparam[1] == 10.0,
			"to->kernelparam[1] is incorrect.");

	gensvm_free_model(from_model);
	gensvm_free_model(to_model);

	return NULL;
}

char *test_copy_model_sigmoid_2()
{
	struct GenModel *from_model = gensvm_init_model();
	struct GenModel *to_model = gensvm_init_model();

	from_model->p = 2.0;
	from_model->lambda = 1.0;
	from_model->epsilon = 1e-8;
	from_model->kappa = 1.0;
	from_model->weight_idx = 2;
	from_model->kerneltype = K_SIGMOID;
	from_model->kernelparam = Calloc(double, 2);
	from_model->kernelparam[0] = 5.0;
	from_model->kernelparam[1] = 10.0;
	to_model->kernelparam = Calloc(double, 2);

	gensvm_copy_model(from_model, to_model);

	mu_assert(to_model->p == 2.0, "to_model->p incorrect.");
	mu_assert(to_model->lambda == 1.0, "to_model->lambda incorrect.");
	mu_assert(to_model->epsilon == 1e-8, "to_model->epsilon incorrect.");
	mu_assert(to_model->kappa == 1.0, "to_model->kappa incorrect.");
	mu_assert(to_model->weight_idx == 2,
		       	"to_model->weight_idx incorrect.");
	mu_assert(to_model->kerneltype == K_SIGMOID,
		       	"to->kerneltype incorrect");
	mu_assert(to_model->kernelparam[0] == 5.0,
			"to->kernelparam[0] is incorrect.");
	mu_assert(to_model->kernelparam[1] == 10.0,
			"to->kernelparam[1] is incorrect.");

	gensvm_free_model(from_model);
	gensvm_free_model(to_model);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_copy_model_linear);
	mu_run_test(test_copy_model_poly_1);
	mu_run_test(test_copy_model_poly_2);
	mu_run_test(test_copy_model_rbf_1);
	mu_run_test(test_copy_model_rbf_2);
	mu_run_test(test_copy_model_sigmoid_1);
	mu_run_test(test_copy_model_sigmoid_2);

	return NULL;
}

RUN_TESTS(all_tests);
