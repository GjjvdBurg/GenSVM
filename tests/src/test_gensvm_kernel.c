/**
 *@file test_gensvm_kernel.c
 *@author G.J.J. van den Burg
 *@date 2016-11-09
 *@brief Unit tests for gensvm_kernel.c
 *
 *@copyright
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
#include "gensvm_kernel.h"

char *test_kernel_copy_kernelparam_to_data_linear()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	model->kerneltype = K_LINEAR;
	data->kerneltype = K_RBF;

	// start test code //
	gensvm_kernel_copy_kernelparam_to_data(model, data);

	mu_assert(data->kerneltype == K_LINEAR, "Incorrect data kerneltype");
	// end test code //

	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_kernel_copy_kernelparam_to_data_rbf()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	model->kerneltype = K_RBF;
	model->gamma = 1.23;

	// start test code //
	gensvm_kernel_copy_kernelparam_to_data(model, data);
	mu_assert(data->kerneltype == K_RBF, "Incorrect data->kerneltype");
	mu_assert(data->gamma == 1.23, "Incorrect data->gamma");

	// end test code //

	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_kernel_copy_kernelparam_to_data_poly()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	model->kerneltype = K_POLY;
	model->gamma = 1.23;
	model->coef = 2.23;
	model->degree = 3.23;

	// start test code //
	gensvm_kernel_copy_kernelparam_to_data(model, data);
	mu_assert(data->kerneltype == K_POLY, "Incorrect data->kerneltype");
	mu_assert(data->gamma == 1.23, "Incorrect data->gamma");
	mu_assert(data->coef == 2.23, "Incorrect data->coef");
	mu_assert(data->degree == 3.23, "Incorrect data->degree");

	// end test code //

	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_kernel_copy_kernelparam_to_data_sigmoid()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	model->kerneltype = K_SIGMOID;
	model->gamma = 1.23;
	model->coef = 2.23;

	// start test code //
	gensvm_kernel_copy_kernelparam_to_data(model, data);
	mu_assert(data->kerneltype == K_SIGMOID, "Incorrect data->kerneltype");
	mu_assert(data->gamma == 1.23, "Incorrect data->gamma");
	mu_assert(data->coef == 2.23, "Incorrect data->coef");

	// end test code //

	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_dot_rbf()
{
	double dot, gamma;
	double *a = Malloc(double, 5);
	double *b = Malloc(double, 5);

	a[0] = 0.5203363837176203;
	a[1] = 0.3860628599460129;
	a[2] = 0.3592536954640216;
	a[3] = 0.6824659760765744;
	a[4] = 0.5390520090020700;

	b[0] = 0.1782643262351465;
	b[1] = 0.0314270210724957;
	b[2] = 0.5887219369641497;
	b[3] = 0.7710042954911620;
	b[4] = 0.8805451245738238;

	// start test code //
	gamma = 1.0;
	dot = gensvm_kernel_dot_rbf(a, b, 5, 1, 1, gamma);
	mu_assert(fabs(dot - 0.657117701533133) < 1e-14, "Incorrect dot (1)");

	gamma = 5.0;
	dot = gensvm_kernel_dot_rbf(a, b, 5, 1, 1, gamma);
	mu_assert(fabs(dot - 0.122522495044048) < 1e-14, "Incorrect dot (2)");

	// end test code //
	free(a);
	free(b);

	return NULL;
}

char *test_dot_poly()
{
	double dot, gamma, coef, degree;
	double *a = Malloc(double, 5);
	double *b = Malloc(double, 5);

	a[0] = 0.5203363837176203;
	a[1] = 0.3860628599460129;
	a[2] = 0.3592536954640216;
	a[3] = 0.6824659760765744;
	a[4] = 0.5390520090020700;

	b[0] = 0.1782643262351465;
	b[1] = 0.0314270210724957;
	b[2] = 0.5887219369641497;
	b[3] = 0.7710042954911620;
	b[4] = 0.8805451245738238;

	// start test code //
	gamma = 1.0;
	coef = 1.0;
	degree = 1.0;
	dot = gensvm_kernel_dot_poly(a, b, 5, 1, 1, gamma, coef, degree);
	mu_assert(fabs(dot - 2.31723456944910) < 1e-14, "Incorrect dot (1)");

	gamma = 1.5;
	coef = 2.5;
	degree = 3.5;
	dot = gensvm_kernel_dot_poly(a, b, 5, 1, 1, gamma, coef, degree);
	mu_assert(fabs(dot - 189.6989652572890179) < 1e-14, "Incorrect dot (2)");

	// end test code //
	free(a);
	free(b);
	return NULL;
}

char *test_dot_sigmoid()
{
	double dot, gamma, coef;
	double *a = Malloc(double, 5);
	double *b = Malloc(double, 5);

	a[0] = 0.5203363837176203;
	a[1] = 0.3860628599460129;
	a[2] = 0.3592536954640216;
	a[3] = 0.6824659760765744;
	a[4] = 0.5390520090020700;

	b[0] = 0.1782643262351465;
	b[1] = 0.0314270210724957;
	b[2] = 0.5887219369641497;
	b[3] = 0.7710042954911620;
	b[4] = 0.8805451245738238;

	// start test code //
	gamma = 1.0;
	coef = 1.0;
	dot = gensvm_kernel_dot_sigmoid(a, b, 5, 1, 1, gamma, coef);
	mu_assert(fabs(dot - 0.9807642810850747) < 1e-14, "Incorrect dot (1)");

	gamma = 1.5;
	coef = 2.5;
	dot = gensvm_kernel_dot_sigmoid(a, b, 5, 1, 1, gamma, coef);
	mu_assert(fabs(dot - 0.9997410009167159) < 1e-14, "Incorrect dot (2)");

	// end test code //
	free(a);
	free(b);

	return NULL;
}


char *test_kernel_preprocess_linear()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();
	model->kerneltype = K_LINEAR;
	data->m = 10;
	data->r = 0;

	gensvm_kernel_preprocess(model, data);

	mu_assert(data->r == 10, "Incorrect GenData::r");

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_kernel_preprocess_kernel()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	model->n = 10;
	model->m = 5;
	data->n = 10;
	data->m = 5;

	model->kerneltype = K_RBF;
	model->gamma = 0.348;
	model->kernel_eigen_cutoff = 5e-3;

	data->Z = Calloc(double, data->n * (data->m + 1));
	data->RAW = data->Z;
	data->Sigma = Calloc(double, 10);

	matrix_set(data->Z, data->n, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data->Z, data->n, data->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data->Z, data->n, data->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data->Z, data->n, data->m+1, 0, 4, 0.8233234072519983);
	matrix_set(data->Z, data->n, data->m+1, 0, 5, 0.2728942849022845);
	matrix_set(data->Z, data->n, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data->Z, data->n, data->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data->Z, data->n, data->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data->Z, data->n, data->m+1, 1, 4, 0.7956168453294307);
	matrix_set(data->Z, data->n, data->m+1, 1, 5, 0.2665038412777220);
	matrix_set(data->Z, data->n, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data->Z, data->n, data->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data->Z, data->n, data->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data->Z, data->n, data->m+1, 2, 4, 0.5777227081256917);
	matrix_set(data->Z, data->n, data->m+1, 2, 5, 0.1126634767453339);
	matrix_set(data->Z, data->n, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data->Z, data->n, data->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data->Z, data->n, data->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data->Z, data->n, data->m+1, 3, 4, 0.4426030703804438);
	matrix_set(data->Z, data->n, data->m+1, 3, 5, 0.9577384085295405);
	matrix_set(data->Z, data->n, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data->Z, data->n, data->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data->Z, data->n, data->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data->Z, data->n, data->m+1, 4, 4, 0.7701104553132680);
	matrix_set(data->Z, data->n, data->m+1, 4, 5, 0.4421972387164297);
	matrix_set(data->Z, data->n, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data->Z, data->n, data->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data->Z, data->n, data->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data->Z, data->n, data->m+1, 5, 4, 0.3267543833513200);
	matrix_set(data->Z, data->n, data->m+1, 5, 5, 0.3993579061297995);
	matrix_set(data->Z, data->n, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data->Z, data->n, data->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data->Z, data->n, data->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data->Z, data->n, data->m+1, 6, 4, 0.3693175185473680);
	matrix_set(data->Z, data->n, data->m+1, 6, 5, 0.9137371110726166);
	matrix_set(data->Z, data->n, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data->Z, data->n, data->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data->Z, data->n, data->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data->Z, data->n, data->m+1, 7, 4, 0.2456426390463990);
	matrix_set(data->Z, data->n, data->m+1, 7, 5, 0.0191704319054138);
	matrix_set(data->Z, data->n, data->m+1, 8, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data->Z, data->n, data->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data->Z, data->n, data->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data->Z, data->n, data->m+1, 8, 4, 0.1102697774064020);
	matrix_set(data->Z, data->n, data->m+1, 8, 5, 0.4826593547737735);
	matrix_set(data->Z, data->n, data->m+1, 9, 0, 1.0000000000000000);
	matrix_set(data->Z, data->n, data->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data->Z, data->n, data->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data->Z, data->n, data->m+1, 9, 3, 0.0884616753393881);
	matrix_set(data->Z, data->n, data->m+1, 9, 4, 0.8659836346403005);
	matrix_set(data->Z, data->n, data->m+1, 9, 5, 0.1028774221216107);

	// start test code //
	gensvm_kernel_preprocess(model, data);
	mu_assert(data->r == 7, "Incorrect data->r");

	double eps = 1e-14;
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 0, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 1)) -
				fabs(2.4632837902141640)) < eps,
			"Incorrect data->Z at 0, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 2)) -
				fabs(-0.3037489220604925)) < eps,
			"Incorrect data->Z at 0, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 3)) -
				fabs(-0.0061287029147240)) < eps,
			"Incorrect data->Z at 0, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 4)) -
				fabs(0.1822712619914593)) < eps,
			"Incorrect data->Z at 0, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 5)) -
				fabs(0.0252737053303148)) < eps,
			"Incorrect data->Z at 0, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 6)) -
				fabs(-0.0078753266252524)) < eps,
			"Incorrect data->Z at 0, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 0, 7)) -
				fabs(-0.0012800124996018)) < eps,
			"Incorrect data->Z at 0, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 1, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 1)) -
				fabs(2.2923640983040641)) < eps,
			"Incorrect data->Z at 1, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 2)) -
				fabs(-0.3048037728463330)) < eps,
			"Incorrect data->Z at 1, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 3)) -
				fabs(-0.2586192720897897)) < eps,
			"Incorrect data->Z at 1, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 4)) -
				fabs(0.1747912247100736)) < eps,
			"Incorrect data->Z at 1, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 5)) -
				fabs(-0.0623497873850738)) < eps,
			"Incorrect data->Z at 1, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 6)) -
				fabs(-0.0199493291395259)) < eps,
			"Incorrect data->Z at 1, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 1, 7)) -
				fabs(0.0068540206892510)) < eps,
			"Incorrect data->Z at 1, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 2, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 1)) -
				fabs(2.4167201742337761)) < eps,
			"Incorrect data->Z at 2, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 2)) -
				fabs(-0.1499385272847361)) < eps,
			"Incorrect data->Z at 2, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 3)) -
				fabs(-0.1781619658696836)) < eps,
			"Incorrect data->Z at 2, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 4)) -
				fabs(-0.2363293887681946)) < eps,
			"Incorrect data->Z at 2, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 5)) -
				fabs(-0.0362117307160720)) < eps,
			"Incorrect data->Z at 2, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 6)) -
				fabs(0.0366137533260933)) < eps,
			"Incorrect data->Z at 2, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 2, 7)) -
				fabs(0.0227049982868101)) < eps,
			"Incorrect data->Z at 2, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 3, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 1)) -
				fabs(2.3339314848494390)) < eps,
			"Incorrect data->Z at 3, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 2)) -
				fabs(0.4067278927345656)) < eps,
			"Incorrect data->Z at 3, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 3)) -
				fabs(0.0198947146890620)) < eps,
			"Incorrect data->Z at 3, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 4)) -
				fabs(0.1187106614859180)) < eps,
			"Incorrect data->Z at 3, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 5)) -
				fabs(0.0734848412140159)) < eps,
			"Incorrect data->Z at 3, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 6)) -
				fabs(-0.0166955533210990)) < eps,
			"Incorrect data->Z at 3, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 3, 7)) -
				fabs(0.0229112619510384)) < eps,
			"Incorrect data->Z at 3, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 4, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 1)) -
				fabs(2.5061509421266424)) < eps,
			"Incorrect data->Z at 4, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 2)) -
				fabs(0.0574469229922174)) < eps,
			"Incorrect data->Z at 4, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 3)) -
				fabs(-0.2858649955147738)) < eps,
			"Incorrect data->Z at 4, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 4)) -
				fabs(-0.0995031375002134)) < eps,
			"Incorrect data->Z at 4, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 5)) -
				fabs(0.0223790101651578)) < eps,
			"Incorrect data->Z at 4, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 6)) -
				fabs(-0.0355571480867735)) < eps,
			"Incorrect data->Z at 4, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 4, 7)) -
				fabs(-0.0219026472149696)) < eps,
			"Incorrect data->Z at 4, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 5, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 1)) -
				fabs(2.4482858151168982)) < eps,
			"Incorrect data->Z at 5, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 2)) -
				fabs(-0.0670998214520230)) < eps,
			"Incorrect data->Z at 5, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 3)) -
				fabs(0.3147064295566219)) < eps,
			"Incorrect data->Z at 5, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 4)) -
				fabs(0.1070535630418465)) < eps,
			"Incorrect data->Z at 5, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 5)) -
				fabs(-0.0052824396955993)) < eps,
			"Incorrect data->Z at 5, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 6)) -
				fabs(0.0614363461130733)) < eps,
			"Incorrect data->Z at 5, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 5, 7)) -
				fabs(-0.0075355247061472)) < eps,
			"Incorrect data->Z at 5, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 6, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 1)) -
				fabs(2.3638928644404329)) < eps,
			"Incorrect data->Z at 6, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 2)) -
				fabs(0.3482541374011597)) < eps,
			"Incorrect data->Z at 6, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 3)) -
				fabs(-0.2422541976251498)) < eps,
			"Incorrect data->Z at 6, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 4)) -
				fabs(0.0251886519764033)) < eps,
			"Incorrect data->Z at 6, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 5)) -
				fabs(-0.0079397861684362)) < eps,
			"Incorrect data->Z at 6, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 6)) -
				fabs(0.0424975213407462)) < eps,
			"Incorrect data->Z at 6, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 6, 7)) -
				fabs(-0.0204279932276333)) < eps,
			"Incorrect data->Z at 6, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 7, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 1)) -
				fabs(2.3607306299135344)) < eps,
			"Incorrect data->Z at 7, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 2)) -
				fabs(-0.0220102589508912)) < eps,
			"Incorrect data->Z at 7, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 3)) -
				fabs(0.3913398731540265)) < eps,
			"Incorrect data->Z at 7, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 4)) -
				fabs(-0.0941469673695446)) < eps,
			"Incorrect data->Z at 7, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 5)) -
				fabs(-0.0477595489009114)) < eps,
			"Incorrect data->Z at 7, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 6)) -
				fabs(-0.0367688438245860)) < eps,
			"Incorrect data->Z at 7, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 7, 7)) -
				fabs(-0.0133498576642393)) < eps,
			"Incorrect data->Z at 7, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 8, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 1)) -
				fabs(2.5023932475093376)) < eps,
			"Incorrect data->Z at 8, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 2)) -
				fabs(0.2929602950386334)) < eps,
			"Incorrect data->Z at 8, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 3)) -
				fabs(0.1371647327912284)) < eps,
			"Incorrect data->Z at 8, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 4)) -
				fabs(-0.0270505649533715)) < eps,
			"Incorrect data->Z at 8, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 5)) -
				fabs(-0.0685258491091892)) < eps,
			"Incorrect data->Z at 8, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 6)) -
				fabs(-0.0213385621647371)) < eps,
			"Incorrect data->Z at 8, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 8, 7)) -
				fabs(0.0121305554343051)) < eps,
			"Incorrect data->Z at 8, 7");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 0)) -
				fabs(1.0000000000000000)) < eps,
			"Incorrect data->Z at 9, 0");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 1)) -
				fabs(2.4579608302226870)) < eps,
			"Incorrect data->Z at 9, 1");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 2)) -
				fabs(-0.2538863526247282)) < eps,
			"Incorrect data->Z at 9, 2");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 3)) -
				fabs(0.0991005899665861)) < eps,
			"Incorrect data->Z at 9, 3");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 4)) -
				fabs(-0.1374781282330359)) < eps,
			"Incorrect data->Z at 9, 4");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 5)) -
				fabs(0.1043628273525485)) < eps,
			"Incorrect data->Z at 9, 5");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 6)) -
				fabs(-0.0024213407513349)) < eps,
			"Incorrect data->Z at 9, 6");
	mu_assert(fabs(fabs(matrix_get(data->Z, data->n, data->r+1, 9, 7)) -
				fabs(0.0007673936590348)) < eps,
			"Incorrect data->Z at 9, 7");

	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 0, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 0, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 0, 1) -
				0.8056271362589000) < eps,
			"Incorrect data->RAW at 0, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 0, 2) -
				0.4874175854113872) < eps,
			"Incorrect data->RAW at 0, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 0, 3) -
				0.4453015882771756) < eps,
			"Incorrect data->RAW at 0, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 0, 4) -
				0.8233234072519983) < eps,
			"Incorrect data->RAW at 0, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 0, 5) -
				0.2728942849022845) < eps,
			"Incorrect data->RAW at 0, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 1, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 1, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 1, 1) -
				0.7940590105180981) < eps,
			"Incorrect data->RAW at 1, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 1, 2) -
				0.1861049005485224) < eps,
			"Incorrect data->RAW at 1, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 1, 3) -
				0.8469394287449229) < eps,
			"Incorrect data->RAW at 1, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 1, 4) -
				0.7956168453294307) < eps,
			"Incorrect data->RAW at 1, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 1, 5) -
				0.2665038412777220) < eps,
			"Incorrect data->RAW at 1, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 2, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 2, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 2, 1) -
				0.0294257611061681) < eps,
			"Incorrect data->RAW at 2, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 2, 2) -
				0.0242717976065267) < eps,
			"Incorrect data->RAW at 2, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 2, 3) -
				0.5039128672814752) < eps,
			"Incorrect data->RAW at 2, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 2, 4) -
				0.5777227081256917) < eps,
			"Incorrect data->RAW at 2, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 2, 5) -
				0.1126634767453339) < eps,
			"Incorrect data->RAW at 2, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 3, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 3, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 3, 1) -
				0.1746563833537603) < eps,
			"Incorrect data->RAW at 3, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 3, 2) -
				0.9135736087631979) < eps,
			"Incorrect data->RAW at 3, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 3, 3) -
				0.5270258081021366) < eps,
			"Incorrect data->RAW at 3, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 3, 4) -
				0.4426030703804438) < eps,
			"Incorrect data->RAW at 3, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 3, 5) -
				0.9577384085295405) < eps,
			"Incorrect data->RAW at 3, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 4, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 4, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 4, 1) -
				0.0022298761599785) < eps,
			"Incorrect data->RAW at 4, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 4, 2) -
				0.3773482059713607) < eps,
			"Incorrect data->RAW at 4, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 4, 3) -
				0.8009654729622842) < eps,
			"Incorrect data->RAW at 4, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 4, 4) -
				0.7701104553132680) < eps,
			"Incorrect data->RAW at 4, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 4, 5) -
				0.4421972387164297) < eps,
			"Incorrect data->RAW at 4, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 5, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 5, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 5, 1) -
				0.6638830667081945) < eps,
			"Incorrect data->RAW at 5, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 5, 2) -
				0.6467607601353914) < eps,
			"Incorrect data->RAW at 5, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 5, 3) -
				0.0434948735457108) < eps,
			"Incorrect data->RAW at 5, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 5, 4) -
				0.3267543833513200) < eps,
			"Incorrect data->RAW at 5, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 5, 5) -
				0.3993579061297995) < eps,
			"Incorrect data->RAW at 5, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 6, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 6, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 6, 1) -
				0.0770493004546461) < eps,
			"Incorrect data->RAW at 6, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 6, 2) -
				0.3699566427075194) < eps,
			"Incorrect data->RAW at 6, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 6, 3) -
				0.7863539761080217) < eps,
			"Incorrect data->RAW at 6, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 6, 4) -
				0.3693175185473680) < eps,
			"Incorrect data->RAW at 6, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 6, 5) -
				0.9137371110726166) < eps,
			"Incorrect data->RAW at 6, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 7, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 7, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 7, 1) -
				0.2685233952731509) < eps,
			"Incorrect data->RAW at 7, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 7, 2) -
				0.8539966432782011) < eps,
			"Incorrect data->RAW at 7, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 7, 3) -
				0.0967159557826836) < eps,
			"Incorrect data->RAW at 7, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 7, 4) -
				0.2456426390463990) < eps,
			"Incorrect data->RAW at 7, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 7, 5) -
				0.0191704319054138) < eps,
			"Incorrect data->RAW at 7, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 8, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 8, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 8, 1) -
				0.1163951898554611) < eps,
			"Incorrect data->RAW at 8, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 8, 2) -
				0.7667861436369238) < eps,
			"Incorrect data->RAW at 8, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 8, 3) -
				0.5031912600213351) < eps,
			"Incorrect data->RAW at 8, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 8, 4) -
				0.1102697774064020) < eps,
			"Incorrect data->RAW at 8, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 8, 5) -
				0.4826593547737735) < eps,
			"Incorrect data->RAW at 8, 5");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 9, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->RAW at 9, 0");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 9, 1) -
				0.2290251898688216) < eps,
			"Incorrect data->RAW at 9, 1");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 9, 2) -
				0.4401981048538806) < eps,
			"Incorrect data->RAW at 9, 2");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 9, 3) -
				0.0884616753393881) < eps,
			"Incorrect data->RAW at 9, 3");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 9, 4) -
				0.8659836346403005) < eps,
			"Incorrect data->RAW at 9, 4");
	mu_assert(fabs(matrix_get(data->RAW, data->n, data->m+1, 9, 5) -
				0.1028774221216107) < eps,
			"Incorrect data->RAW at 9, 5");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}


char *test_kernel_postprocess_linear()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	model->kerneltype = K_LINEAR;
	test->m = 10;
	test->r = 0;

	// start test code //
	gensvm_kernel_postprocess(model, train, test);

	mu_assert(test->r == 10, "Incorrect test->r");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(train);
	gensvm_free_data(test);

	return NULL;
}

char *test_kernel_postprocess_kernel()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	train->n = 10;
	train->m = 5;
	train->r = train->m;
	train->RAW = Calloc(double, train->n * (train->m + 1));
	train->Sigma = Calloc(double, train->m);

	test->n = 8;
	test->m = 5;
	test->RAW = Calloc(double, test->n * (test->m + 1));

	model->kerneltype = K_RBF;
	model->gamma = 1.132;

	// start test code //

	matrix_set(train->RAW, train->n, train->m+1, 0, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 0, 1, 0.8056271362589000);
	matrix_set(train->RAW, train->n, train->m+1, 0, 2, 0.4874175854113872);
	matrix_set(train->RAW, train->n, train->m+1, 0, 3, 0.4453015882771756);
	matrix_set(train->RAW, train->n, train->m+1, 0, 4, 0.8233234072519983);
	matrix_set(train->RAW, train->n, train->m+1, 0, 5, 0.2728942849022845);
	matrix_set(train->RAW, train->n, train->m+1, 1, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 1, 1, 0.7940590105180981);
	matrix_set(train->RAW, train->n, train->m+1, 1, 2, 0.1861049005485224);
	matrix_set(train->RAW, train->n, train->m+1, 1, 3, 0.8469394287449229);
	matrix_set(train->RAW, train->n, train->m+1, 1, 4, 0.7956168453294307);
	matrix_set(train->RAW, train->n, train->m+1, 1, 5, 0.2665038412777220);
	matrix_set(train->RAW, train->n, train->m+1, 2, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 2, 1, 0.0294257611061681);
	matrix_set(train->RAW, train->n, train->m+1, 2, 2, 0.0242717976065267);
	matrix_set(train->RAW, train->n, train->m+1, 2, 3, 0.5039128672814752);
	matrix_set(train->RAW, train->n, train->m+1, 2, 4, 0.5777227081256917);
	matrix_set(train->RAW, train->n, train->m+1, 2, 5, 0.1126634767453339);
	matrix_set(train->RAW, train->n, train->m+1, 3, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 3, 1, 0.1746563833537603);
	matrix_set(train->RAW, train->n, train->m+1, 3, 2, 0.9135736087631979);
	matrix_set(train->RAW, train->n, train->m+1, 3, 3, 0.5270258081021366);
	matrix_set(train->RAW, train->n, train->m+1, 3, 4, 0.4426030703804438);
	matrix_set(train->RAW, train->n, train->m+1, 3, 5, 0.9577384085295405);
	matrix_set(train->RAW, train->n, train->m+1, 4, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 4, 1, 0.0022298761599785);
	matrix_set(train->RAW, train->n, train->m+1, 4, 2, 0.3773482059713607);
	matrix_set(train->RAW, train->n, train->m+1, 4, 3, 0.8009654729622842);
	matrix_set(train->RAW, train->n, train->m+1, 4, 4, 0.7701104553132680);
	matrix_set(train->RAW, train->n, train->m+1, 4, 5, 0.4421972387164297);
	matrix_set(train->RAW, train->n, train->m+1, 5, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 5, 1, 0.6638830667081945);
	matrix_set(train->RAW, train->n, train->m+1, 5, 2, 0.6467607601353914);
	matrix_set(train->RAW, train->n, train->m+1, 5, 3, 0.0434948735457108);
	matrix_set(train->RAW, train->n, train->m+1, 5, 4, 0.3267543833513200);
	matrix_set(train->RAW, train->n, train->m+1, 5, 5, 0.3993579061297995);
	matrix_set(train->RAW, train->n, train->m+1, 6, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 6, 1, 0.0770493004546461);
	matrix_set(train->RAW, train->n, train->m+1, 6, 2, 0.3699566427075194);
	matrix_set(train->RAW, train->n, train->m+1, 6, 3, 0.7863539761080217);
	matrix_set(train->RAW, train->n, train->m+1, 6, 4, 0.3693175185473680);
	matrix_set(train->RAW, train->n, train->m+1, 6, 5, 0.9137371110726166);
	matrix_set(train->RAW, train->n, train->m+1, 7, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 7, 1, 0.2685233952731509);
	matrix_set(train->RAW, train->n, train->m+1, 7, 2, 0.8539966432782011);
	matrix_set(train->RAW, train->n, train->m+1, 7, 3, 0.0967159557826836);
	matrix_set(train->RAW, train->n, train->m+1, 7, 4, 0.2456426390463990);
	matrix_set(train->RAW, train->n, train->m+1, 7, 5, 0.0191704319054138);
	matrix_set(train->RAW, train->n, train->m+1, 8, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 8, 1, 0.1163951898554611);
	matrix_set(train->RAW, train->n, train->m+1, 8, 2, 0.7667861436369238);
	matrix_set(train->RAW, train->n, train->m+1, 8, 3, 0.5031912600213351);
	matrix_set(train->RAW, train->n, train->m+1, 8, 4, 0.1102697774064020);
	matrix_set(train->RAW, train->n, train->m+1, 8, 5, 0.4826593547737735);
	matrix_set(train->RAW, train->n, train->m+1, 9, 0, 1.0000000000000000);
	matrix_set(train->RAW, train->n, train->m+1, 9, 1, 0.2290251898688216);
	matrix_set(train->RAW, train->n, train->m+1, 9, 2, 0.4401981048538806);
	matrix_set(train->RAW, train->n, train->m+1, 9, 3, 0.0884616753393881);
	matrix_set(train->RAW, train->n, train->m+1, 9, 4, 0.8659836346403005);
	matrix_set(train->RAW, train->n, train->m+1, 9, 5, 0.1028774221216107);
	train->Z = train->RAW;

	matrix_set(test->RAW, test->n, test->m+1, 0, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 0, 1, 0.4221779040953508);
	matrix_set(test->RAW, test->n, test->m+1, 0, 2, 0.2491297759535341);
	matrix_set(test->RAW, test->n, test->m+1, 0, 3, 0.3534695686769205);
	matrix_set(test->RAW, test->n, test->m+1, 0, 4, 0.1843977014611545);
	matrix_set(test->RAW, test->n, test->m+1, 0, 5, 0.2520805490786355);
	matrix_set(test->RAW, test->n, test->m+1, 1, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 1, 1, 0.1657896146622650);
	matrix_set(test->RAW, test->n, test->m+1, 1, 2, 0.6189568847484389);
	matrix_set(test->RAW, test->n, test->m+1, 1, 3, 0.0966039750517987);
	matrix_set(test->RAW, test->n, test->m+1, 1, 4, 0.2271271152279254);
	matrix_set(test->RAW, test->n, test->m+1, 1, 5, 0.6868351025871483);
	matrix_set(test->RAW, test->n, test->m+1, 2, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 2, 1, 0.2472125660924395);
	matrix_set(test->RAW, test->n, test->m+1, 2, 2, 0.5090230605656835);
	matrix_set(test->RAW, test->n, test->m+1, 2, 3, 0.7277927065471604);
	matrix_set(test->RAW, test->n, test->m+1, 2, 4, 0.6985271697744179);
	matrix_set(test->RAW, test->n, test->m+1, 2, 5, 0.6358202467851303);
	matrix_set(test->RAW, test->n, test->m+1, 3, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 3, 1, 0.7987328321610475);
	matrix_set(test->RAW, test->n, test->m+1, 3, 2, 0.9000985899630121);
	matrix_set(test->RAW, test->n, test->m+1, 3, 3, 0.6893982492936992);
	matrix_set(test->RAW, test->n, test->m+1, 3, 4, 0.3663483426181376);
	matrix_set(test->RAW, test->n, test->m+1, 3, 5, 0.6240335158701935);
	matrix_set(test->RAW, test->n, test->m+1, 4, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 4, 1, 0.9763647809741268);
	matrix_set(test->RAW, test->n, test->m+1, 4, 2, 0.1140159498084029);
	matrix_set(test->RAW, test->n, test->m+1, 4, 3, 0.8917076715220298);
	matrix_set(test->RAW, test->n, test->m+1, 4, 4, 0.0724688230975385);
	matrix_set(test->RAW, test->n, test->m+1, 4, 5, 0.3762619556957781);
	matrix_set(test->RAW, test->n, test->m+1, 5, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 5, 1, 0.5095430570208566);
	matrix_set(test->RAW, test->n, test->m+1, 5, 2, 0.3920414766506218);
	matrix_set(test->RAW, test->n, test->m+1, 5, 3, 0.7448958325189663);
	matrix_set(test->RAW, test->n, test->m+1, 5, 4, 0.5028223670359105);
	matrix_set(test->RAW, test->n, test->m+1, 5, 5, 0.0659571469966666);
	matrix_set(test->RAW, test->n, test->m+1, 6, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 6, 1, 0.1190008864254177);
	matrix_set(test->RAW, test->n, test->m+1, 6, 2, 0.7735127328107272);
	matrix_set(test->RAW, test->n, test->m+1, 6, 3, 0.7292497207505059);
	matrix_set(test->RAW, test->n, test->m+1, 6, 4, 0.3385406475521640);
	matrix_set(test->RAW, test->n, test->m+1, 6, 5, 0.7274062216238782);
	matrix_set(test->RAW, test->n, test->m+1, 7, 0, 1.0000000000000000);
	matrix_set(test->RAW, test->n, test->m+1, 7, 1, 0.5211738041619068);
	matrix_set(test->RAW, test->n, test->m+1, 7, 2, 0.0437246682449493);
	matrix_set(test->RAW, test->n, test->m+1, 7, 3, 0.0388180555183835);
	matrix_set(test->RAW, test->n, test->m+1, 7, 4, 0.6866094956337303);
	matrix_set(test->RAW, test->n, test->m+1, 7, 5, 0.7166516653360046);
	test->Z = test->RAW;

	matrix_set(train->Sigma, train->m, 1, 0, 0, 0.9489494575279929);
	matrix_set(train->Sigma, train->m, 1, 1, 0, 0.5162086311506285);
	matrix_set(train->Sigma, train->m, 1, 2, 0, 0.6279732910944019);
	matrix_set(train->Sigma, train->m, 1, 3, 0, 0.4247685972225844);
	matrix_set(train->Sigma, train->m, 1, 4, 0, 0.9197312157223199);

	gensvm_kernel_postprocess(model, train, test);

	// start test code //
	double eps = 1e-14;
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 0, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 1) -
				1.7919375085966578) < eps,
			"Incorrect test->Z at 0, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 2) -
				9.3749010505955841) < eps,
			"Incorrect test->Z at 0, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 3) -
				5.4188812347396071) < eps,
			"Incorrect test->Z at 0, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 4) -
				14.1749662709433490) < eps,
			"Incorrect test->Z at 0, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 5) -
				2.0967148208956785) < eps,
			"Incorrect test->Z at 0, 5");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 1, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 1) -
				1.4492639390367632) < eps,
			"Incorrect test->Z at 1, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 2) -
				10.3598279425411341) < eps,
			"Incorrect test->Z at 1, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 3) -
				4.9086179135381158) < eps,
			"Incorrect test->Z at 1, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 4) -
				11.5763085861210531) < eps,
			"Incorrect test->Z at 1, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 5) -
				2.4474439570061368) < eps,
			"Incorrect test->Z at 1, 5");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 2, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 1) -
				1.7125210555223442) < eps,
			"Incorrect test->Z at 2, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 2) -
				10.1267642728449463) < eps,
			"Incorrect test->Z at 2, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 3) -
				7.4711312071081304) < eps,
			"Incorrect test->Z at 2, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 4) -
				16.7954195638576316) < eps,
			"Incorrect test->Z at 2, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 5) -
				2.9867291329866701) < eps,
			"Incorrect test->Z at 2, 5");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 3, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 1) -
				1.5599347974839732) < eps,
			"Incorrect test->Z at 3, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 2) -
				8.1287657651835250) < eps,
			"Incorrect test->Z at 3, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 3) -
				4.4627669627170601) < eps,
			"Incorrect test->Z at 3, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 4) -
				10.5382257909839403) < eps,
			"Incorrect test->Z at 3, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 5) -
				2.0425829022313486) < eps,
			"Incorrect test->Z at 3, 5");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 4, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 1) -
				1.1070615251210820) < eps,
			"Incorrect test->Z at 4, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 2) -
				3.8591733667970094) < eps,
			"Incorrect test->Z at 4, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 3) -
				3.1420027555179613) < eps,
			"Incorrect test->Z at 4, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 4) -
				7.3657520580430189) < eps,
			"Incorrect test->Z at 4, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 5) -
				1.0570891668650841) < eps,
			"Incorrect test->Z at 4, 5");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 5, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 5, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 5, 1) -
				2.0048972078646314) < eps,
			"Incorrect test->Z at 5, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 5, 2) -
				8.5523796391619555) < eps,
			"Incorrect test->Z at 5, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 5, 3) -
				6.2297195520069808) < eps,
			"Incorrect test->Z at 5, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 5, 4) -
				16.1775074795210649) < eps,
			"Incorrect test->Z at 5, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 5, 5) -
				2.0441111634728997) < eps,
			"Incorrect test->Z at 5, 5");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 6, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 6, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 6, 1) -
				1.2742361109643845) < eps,
			"Incorrect test->Z at 6, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 6, 2) -
				10.4770471122282736) < eps,
			"Incorrect test->Z at 6, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 6, 3) -
				6.4647232181010859) < eps,
			"Incorrect test->Z at 6, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 6, 4) -
				12.6115089203184212) < eps,
			"Incorrect test->Z at 6, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 6, 5) -
				3.0535665607782896) < eps,
			"Incorrect test->Z at 6, 5");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 7, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 7, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 7, 1) -
				1.3901612934200134) < eps,
			"Incorrect test->Z at 7, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 7, 2) -
				6.1897209583367783) < eps,
			"Incorrect test->Z at 7, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 7, 3) -
				3.8589726502511379) < eps,
			"Incorrect test->Z at 7, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 7, 4) -
				11.0809129571265199) < eps,
			"Incorrect test->Z at 7, 4");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 7, 5) -
				1.5298881436115606) < eps,
			"Incorrect test->Z at 7, 5");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(train);
	gensvm_free_data(test);

	return NULL;
}

char *test_kernel_compute_rbf()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	// setup //
	data->n = 10;
	data->m = 3;
	data->RAW = Calloc(double, data->n * (data->m + 1));

	model->n = 10;
	model->m = 3;
	model->kerneltype = K_RBF;
	model->gamma = 0.348;

	double *K = Calloc(double, data->n * data->n);

	matrix_set(data->RAW, data->n, data->m+1, 0, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data->RAW, data->n, data->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data->RAW, data->n, data->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data->RAW, data->n, data->m+1, 1, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data->RAW, data->n, data->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data->RAW, data->n, data->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data->RAW, data->n, data->m+1, 2, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data->RAW, data->n, data->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data->RAW, data->n, data->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data->RAW, data->n, data->m+1, 3, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data->RAW, data->n, data->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data->RAW, data->n, data->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data->RAW, data->n, data->m+1, 4, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data->RAW, data->n, data->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data->RAW, data->n, data->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data->RAW, data->n, data->m+1, 5, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data->RAW, data->n, data->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data->RAW, data->n, data->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data->RAW, data->n, data->m+1, 6, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data->RAW, data->n, data->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data->RAW, data->n, data->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data->RAW, data->n, data->m+1, 7, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data->RAW, data->n, data->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data->RAW, data->n, data->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data->RAW, data->n, data->m+1, 8, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data->RAW, data->n, data->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data->RAW, data->n, data->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data->RAW, data->n, data->m+1, 9, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data->RAW, data->n, data->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data->RAW, data->n, data->m+1, 9, 3, 0.0884616753393881);

	// end setup //

	// start test code //
	double eps = 1e-14;
	gensvm_kernel_compute(model, data, K);

	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 0) -
				1.0000000000000000) < eps,
			"Incorrect K at 0, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 1) -
				0.9159640457437349) < eps,
			"Incorrect K at 0, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 2) -
				0.7516316505770226) < eps,
			"Incorrect K at 0, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 3) -
				0.8154041985130495) < eps,
			"Incorrect K at 0, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 4) -
				0.7612049732149070) < eps,
			"Incorrect K at 0, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 5) -
				0.9305199788157297) < eps,
			"Incorrect K at 0, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 6) -
				0.7945280927801784) < eps,
			"Incorrect K at 0, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 7) -
				0.8274209897573236) < eps,
			"Incorrect K at 0, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 8) -
				0.8239539297657018) < eps,
			"Incorrect K at 0, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 9) -
				0.8514726457528513) < eps,
			"Incorrect K at 0, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 0) -
				0.9159640457437349) < eps,
			"Incorrect K at 1, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 1) -
				1.0000000000000000) < eps,
			"Incorrect K at 1, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 2) -
				0.7760600552154528) < eps,
			"Incorrect K at 1, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 3) -
				0.7023700752771812) < eps,
			"Incorrect K at 1, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 4) -
				0.7932181642748459) < eps,
			"Incorrect K at 1, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 5) -
				0.7375760206683935) < eps,
			"Incorrect K at 1, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 6) -
				0.8253497074841124) < eps,
			"Incorrect K at 1, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 7) -
				0.6394060322576254) < eps,
			"Incorrect K at 1, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 8) -
				0.7274031225620441) < eps,
			"Incorrect K at 1, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 9) -
				0.7162170274961248) < eps,
			"Incorrect K at 1, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 0) -
				0.7516316505770226) < eps,
			"Incorrect K at 2, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 1) -
				0.7760600552154528) < eps,
			"Incorrect K at 2, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 2) -
				1.0000000000000000) < eps,
			"Incorrect K at 2, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 3) -
				0.7537124265723788) < eps,
			"Incorrect K at 2, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 4) -
				0.9283488505800269) < eps,
			"Incorrect K at 2, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 5) -
				0.7056043340833220) < eps,
			"Incorrect K at 2, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 6) -
				0.9322674240595221) < eps,
			"Incorrect K at 2, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 7) -
				0.7282038382529162) < eps,
			"Incorrect K at 2, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 8) -
				0.8232508450896556) < eps,
			"Incorrect K at 2, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 9) -
				0.8744753760805232) < eps,
			"Incorrect K at 2, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 0) -
				0.8154041985130495) < eps,
			"Incorrect K at 3, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 1) -
				0.7023700752771812) < eps,
			"Incorrect K at 3, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 2) -
				0.7537124265723788) < eps,
			"Incorrect K at 3, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 3) -
				1.0000000000000000) < eps,
			"Incorrect K at 3, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 4) -
				0.8723850250880842) < eps,
			"Incorrect K at 3, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 5) -
				0.8274320954451651) < eps,
			"Incorrect K at 3, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 6) -
				0.8784822458273475) < eps,
			"Incorrect K at 3, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 7) -
				0.9335699405227060) < eps,
			"Incorrect K at 3, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 8) -
				0.9911621491714698) < eps,
			"Incorrect K at 3, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 9) -
				0.8642062408184789) < eps,
			"Incorrect K at 3, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 0) -
				0.7612049732149070) < eps,
			"Incorrect K at 4, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 1) -
				0.7932181642748459) < eps,
			"Incorrect K at 4, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 2) -
				0.9283488505800269) < eps,
			"Incorrect K at 4, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 3) -
				0.8723850250880842) < eps,
			"Incorrect K at 4, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 4) -
				1.0000000000000000) < eps,
			"Incorrect K at 4, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 5) -
				0.6857259898857532) < eps,
			"Incorrect K at 4, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 6) -
				0.9979606873291865) < eps,
			"Incorrect K at 4, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 7) -
				0.7585568540994299) < eps,
			"Incorrect K at 4, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 8) -
				0.9156042422655802) < eps,
			"Incorrect K at 4, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 9) -
				0.8220610108488617) < eps,
			"Incorrect K at 4, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 0) -
				0.9305199788157297) < eps,
			"Incorrect K at 5, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 1) -
				0.7375760206683935) < eps,
			"Incorrect K at 5, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 2) -
				0.7056043340833220) < eps,
			"Incorrect K at 5, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 3) -
				0.8274320954451651) < eps,
			"Incorrect K at 5, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 4) -
				0.6857259898857532) < eps,
			"Incorrect K at 5, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 5) -
				1.0000000000000000) < eps,
			"Incorrect K at 5, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 6) -
				0.7128058155310608) < eps,
			"Incorrect K at 5, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 7) -
				0.9320891767185483) < eps,
			"Incorrect K at 5, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 8) -
				0.8328818339901494) < eps,
			"Incorrect K at 5, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 9) -
				0.9218622704537242) < eps,
			"Incorrect K at 5, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 0) -
				0.7945280927801784) < eps,
			"Incorrect K at 6, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 1) -
				0.8253497074841124) < eps,
			"Incorrect K at 6, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 2) -
				0.9322674240595221) < eps,
			"Incorrect K at 6, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 3) -
				0.8784822458273475) < eps,
			"Incorrect K at 6, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 4) -
				0.9979606873291865) < eps,
			"Incorrect K at 6, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 5) -
				0.7128058155310608) < eps,
			"Incorrect K at 6, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 6) -
				1.0000000000000000) < eps,
			"Incorrect K at 6, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 7) -
				0.7712042211054879) < eps,
			"Incorrect K at 6, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 8) -
				0.9201278964602475) < eps,
			"Incorrect K at 6, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 9) -
				0.8358974835836517) < eps,
			"Incorrect K at 6, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 0) -
				0.8274209897573236) < eps,
			"Incorrect K at 7, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 1) -
				0.6394060322576254) < eps,
			"Incorrect K at 7, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 2) -
				0.7282038382529162) < eps,
			"Incorrect K at 7, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 3) -
				0.9335699405227060) < eps,
			"Incorrect K at 7, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 4) -
				0.7585568540994299) < eps,
			"Incorrect K at 7, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 5) -
				0.9320891767185483) < eps,
			"Incorrect K at 7, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 6) -
				0.7712042211054879) < eps,
			"Incorrect K at 7, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 7) -
				1.0000000000000000) < eps,
			"Incorrect K at 7, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 8) -
				0.9340756478568407) < eps,
			"Incorrect K at 7, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 9) -
				0.9416191361969405) < eps,
			"Incorrect K at 7, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 0) -
				0.8239539297657018) < eps,
			"Incorrect K at 8, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 1) -
				0.7274031225620441) < eps,
			"Incorrect K at 8, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 2) -
				0.8232508450896556) < eps,
			"Incorrect K at 8, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 3) -
				0.9911621491714698) < eps,
			"Incorrect K at 8, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 4) -
				0.9156042422655802) < eps,
			"Incorrect K at 8, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 5) -
				0.8328818339901494) < eps,
			"Incorrect K at 8, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 6) -
				0.9201278964602475) < eps,
			"Incorrect K at 8, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 7) -
				0.9340756478568407) < eps,
			"Incorrect K at 8, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 8) -
				1.0000000000000000) < eps,
			"Incorrect K at 8, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 9) -
				0.9035820400773296) < eps,
			"Incorrect K at 8, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 0) -
				0.8514726457528513) < eps,
			"Incorrect K at 9, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 1) -
				0.7162170274961248) < eps,
			"Incorrect K at 9, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 2) -
				0.8744753760805232) < eps,
			"Incorrect K at 9, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 3) -
				0.8642062408184789) < eps,
			"Incorrect K at 9, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 4) -
				0.8220610108488617) < eps,
			"Incorrect K at 9, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 5) -
				0.9218622704537242) < eps,
			"Incorrect K at 9, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 6) -
				0.8358974835836517) < eps,
			"Incorrect K at 9, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 7) -
				0.9416191361969405) < eps,
			"Incorrect K at 9, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 8) -
				0.9035820400773296) < eps,
			"Incorrect K at 9, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 9) -
				1.0000000000000000) < eps,
			"Incorrect K at 9, 9");

	// end test code //

	free(K);
	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_kernel_compute_poly()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	// setup //
	data->n = 10;
	data->m = 3;
	data->RAW = Calloc(double, data->n * (data->m + 1));

	model->n = 10;
	model->m = 3;
	model->kerneltype = K_POLY;
	model->gamma = 1.5;
	model->coef = 3.0;
	model->degree = 1.78;

	double *K = Calloc(double, data->n * data->n);

	matrix_set(data->RAW, data->n, data->m+1, 0, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data->RAW, data->n, data->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data->RAW, data->n, data->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data->RAW, data->n, data->m+1, 1, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data->RAW, data->n, data->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data->RAW, data->n, data->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data->RAW, data->n, data->m+1, 2, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data->RAW, data->n, data->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data->RAW, data->n, data->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data->RAW, data->n, data->m+1, 3, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data->RAW, data->n, data->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data->RAW, data->n, data->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data->RAW, data->n, data->m+1, 4, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data->RAW, data->n, data->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data->RAW, data->n, data->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data->RAW, data->n, data->m+1, 5, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data->RAW, data->n, data->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data->RAW, data->n, data->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data->RAW, data->n, data->m+1, 6, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data->RAW, data->n, data->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data->RAW, data->n, data->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data->RAW, data->n, data->m+1, 7, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data->RAW, data->n, data->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data->RAW, data->n, data->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data->RAW, data->n, data->m+1, 8, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data->RAW, data->n, data->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data->RAW, data->n, data->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data->RAW, data->n, data->m+1, 9, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data->RAW, data->n, data->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data->RAW, data->n, data->m+1, 9, 3, 0.0884616753393881);

	// end setup //

	// start test code //
	double eps = 1e-14;
	gensvm_kernel_compute(model, data, K);

	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 0) -
				15.2859595967658617) < eps,
			"Incorrect K at 0, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 1) -
				15.4864409009677200) < eps,
			"Incorrect K at 0, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 2) -
				8.7847919771675169) < eps,
			"Incorrect K at 0, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 3) -
				13.0338783338890352) < eps,
			"Incorrect K at 0, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 4) -
				10.8336465195648390) < eps,
			"Incorrect K at 0, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 5) -
				13.4376986903046944) < eps,
			"Incorrect K at 0, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 6) -
				11.2170859716846412) < eps,
			"Incorrect K at 0, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 7) -
				11.8649726239060982) < eps,
			"Incorrect K at 0, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 8) -
				11.9910642907247134) < eps,
			"Incorrect K at 0, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 9) -
				10.0579319882111839) < eps,
			"Incorrect K at 0, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 0) -
				15.4864409009677200) < eps,
			"Incorrect K at 1, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 1) -
				18.0085663057827681) < eps,
			"Incorrect K at 1, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 2) -
				10.1772186353572813) < eps,
			"Incorrect K at 1, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 3) -
				12.4990693290747981) < eps,
			"Incorrect K at 1, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 4) -
				12.4611344803924595) < eps,
			"Incorrect K at 1, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 5) -
				11.9338476573364556) < eps,
			"Incorrect K at 1, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 6) -
				12.8317520536684597) < eps,
			"Incorrect K at 1, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 7) -
				10.1728085503000116) < eps,
			"Incorrect K at 1, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 8) -
				11.7519359337010130) < eps,
			"Incorrect K at 1, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 9) -
				9.3372389753548415) < eps,
			"Incorrect K at 1, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 0) -
				8.7847919771675169) < eps,
			"Incorrect K at 2, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 1) -
				10.1772186353572813) < eps,
			"Incorrect K at 2, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 2) -
				8.7533569584420832) < eps,
			"Incorrect K at 2, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 3) -
				9.0141355927477953) < eps,
			"Incorrect K at 2, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 4) -
				9.8706020767326290) < eps,
			"Incorrect K at 2, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 5) -
				7.4311529847899962) < eps,
			"Incorrect K at 2, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 6) -
				9.8317460227342295) < eps,
			"Incorrect K at 2, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 7) -
				7.5616330506213130) < eps,
			"Incorrect K at 2, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 8) -
				8.8935156801233308) < eps,
			"Incorrect K at 2, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 9) -
				7.4623650201501697) < eps,
			"Incorrect K at 2, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 0) -
				13.0338783338890352) < eps,
			"Incorrect K at 3, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 1) -
				12.4990693290747981) < eps,
			"Incorrect K at 3, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 2) -
				9.0141355927477953) < eps,
			"Incorrect K at 3, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 3) -
				15.8010293983089412) < eps,
			"Incorrect K at 3, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 4) -
				12.5976572300545939) < eps,
			"Incorrect K at 3, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 5) -
				12.2952603928598840) < eps,
			"Incorrect K at 3, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 6) -
				12.5864308488164252) < eps,
			"Incorrect K at 3, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 7) -
				13.5095269679848702) < eps,
			"Incorrect K at 3, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 8) -
				14.4248904770779447) < eps,
			"Incorrect K at 3, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 9) -
				10.4303067289663023) < eps,
			"Incorrect K at 3, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 0) -
				10.8336465195648390) < eps,
			"Incorrect K at 4, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 1) -
				12.4611344803924595) < eps,
			"Incorrect K at 4, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 2) -
				9.8706020767326290) < eps,
			"Incorrect K at 4, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 3) -
				12.5976572300545939) < eps,
			"Incorrect K at 4, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 4) -
				12.7332028710821596) < eps,
			"Incorrect K at 4, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 5) -
				8.9267358166744923) < eps,
			"Incorrect K at 4, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 6) -
				12.6168065443066961) < eps,
			"Incorrect K at 4, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 7) -
				9.7796024413467819) < eps,
			"Incorrect K at 4, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 8) -
				11.9994372472037618) < eps,
			"Incorrect K at 4, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 9) -
				8.6300134557803005) < eps,
			"Incorrect K at 4, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 0) -
				13.4376986903046944) < eps,
			"Incorrect K at 5, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 1) -
				11.9338476573364556) < eps,
			"Incorrect K at 5, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 2) -
				7.4311529847899962) < eps,
			"Incorrect K at 5, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 3) -
				12.2952603928598840) < eps,
			"Incorrect K at 5, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 4) -
				8.9267358166744923) < eps,
			"Incorrect K at 5, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 5) -
				13.3667509754465605) < eps,
			"Incorrect K at 5, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 6) -
				9.2374362113086974) < eps,
			"Incorrect K at 5, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 7) -
				12.3359255360486753) < eps,
			"Incorrect K at 5, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 8) -
				11.2365407617438784) < eps,
			"Incorrect K at 5, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 9) -
				10.0736132325944858) < eps,
			"Incorrect K at 5, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 0) -
				11.2170859716846412) < eps,
			"Incorrect K at 6, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 1) -
				12.8317520536684597) < eps,
			"Incorrect K at 6, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 2) -
				9.8317460227342295) < eps,
			"Incorrect K at 6, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 3) -
				12.5864308488164252) < eps,
			"Incorrect K at 6, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 4) -
				12.6168065443066961) < eps,
			"Incorrect K at 6, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 5) -
				9.2374362113086974) < eps,
			"Incorrect K at 6, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 6) -
				12.5482922903136540) < eps,
			"Incorrect K at 6, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 7) -
				9.8694581744218173) < eps,
			"Incorrect K at 6, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 8) -
				11.9652582335280364) < eps,
			"Incorrect K at 6, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 9) -
				8.7166346335015064) < eps,
			"Incorrect K at 6, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 0) -
				11.8649726239060982) < eps,
			"Incorrect K at 7, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 1) -
				10.1728085503000116) < eps,
			"Incorrect K at 7, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 2) -
				7.5616330506213130) < eps,
			"Incorrect K at 7, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 3) -
				13.5095269679848702) < eps,
			"Incorrect K at 7, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 4) -
				9.7796024413467819) < eps,
			"Incorrect K at 7, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 5) -
				12.3359255360486753) < eps,
			"Incorrect K at 7, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 6) -
				9.8694581744218173) < eps,
			"Incorrect K at 7, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 7) -
				12.9524293684505221) < eps,
			"Incorrect K at 7, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 8) -
				12.3355075975696309) < eps,
			"Incorrect K at 7, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 9) -
				10.1132148127652659) < eps,
			"Incorrect K at 7, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 0) -
				11.9910642907247134) < eps,
			"Incorrect K at 8, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 1) -
				11.7519359337010130) < eps,
			"Incorrect K at 8, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 2) -
				8.8935156801233308) < eps,
			"Incorrect K at 8, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 3) -
				14.4248904770779447) < eps,
			"Incorrect K at 8, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 4) -
				11.9994372472037618) < eps,
			"Incorrect K at 8, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 5) -
				11.2365407617438784) < eps,
			"Incorrect K at 8, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 6) -
				11.9652582335280364) < eps,
			"Incorrect K at 8, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 7) -
				12.3355075975696309) < eps,
			"Incorrect K at 8, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 8) -
				13.3150520649955109) < eps,
			"Incorrect K at 8, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 9) -
				9.8405270375666198) < eps,
			"Incorrect K at 8, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 0) -
				10.0579319882111839) < eps,
			"Incorrect K at 9, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 1) -
				9.3372389753548415) < eps,
			"Incorrect K at 9, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 2) -
				7.4623650201501697) < eps,
			"Incorrect K at 9, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 3) -
				10.4303067289663023) < eps,
			"Incorrect K at 9, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 4) -
				8.6300134557803005) < eps,
			"Incorrect K at 9, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 5) -
				10.0736132325944858) < eps,
			"Incorrect K at 9, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 6) -
				8.7166346335015064) < eps,
			"Incorrect K at 9, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 7) -
				10.1132148127652659) < eps,
			"Incorrect K at 9, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 8) -
				9.8405270375666198) < eps,
			"Incorrect K at 9, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 9) -
				8.7441654683119356) < eps,
			"Incorrect K at 9, 9");
	// end test code //

	free(K);
	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_kernel_compute_sigmoid()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	// setup //
	data->n = 10;
	data->m = 3;
	data->RAW = Calloc(double, data->n * (data->m + 1));

	model->n = 10;
	model->m = 3;
	model->kerneltype = K_SIGMOID;
	model->gamma = 1.23;
	model->coef = 1.6;

	double *K = Calloc(double, data->n * data->n);

	matrix_set(data->RAW, data->n, data->m+1, 0, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data->RAW, data->n, data->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data->RAW, data->n, data->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data->RAW, data->n, data->m+1, 1, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data->RAW, data->n, data->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data->RAW, data->n, data->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data->RAW, data->n, data->m+1, 2, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data->RAW, data->n, data->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data->RAW, data->n, data->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data->RAW, data->n, data->m+1, 3, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data->RAW, data->n, data->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data->RAW, data->n, data->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data->RAW, data->n, data->m+1, 4, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data->RAW, data->n, data->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data->RAW, data->n, data->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data->RAW, data->n, data->m+1, 5, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data->RAW, data->n, data->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data->RAW, data->n, data->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data->RAW, data->n, data->m+1, 6, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data->RAW, data->n, data->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data->RAW, data->n, data->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data->RAW, data->n, data->m+1, 7, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data->RAW, data->n, data->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data->RAW, data->n, data->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data->RAW, data->n, data->m+1, 8, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data->RAW, data->n, data->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data->RAW, data->n, data->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data->RAW, data->n, data->m+1, 9, 0, 1.0);
	matrix_set(data->RAW, data->n, data->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data->RAW, data->n, data->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data->RAW, data->n, data->m+1, 9, 3, 0.0884616753393881);

	// end setup //

	// start test code //
	double eps = 1e-14;
	gensvm_kernel_compute(model, data, K);

	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 0) -
				0.9943637704548856) < eps,
			"Incorrect K at 0, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 1) -
				0.9946686121177555) < eps,
			"Incorrect K at 0, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 2) -
				0.9578937155783831) < eps,
			"Incorrect K at 0, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 3) -
				0.9892315762620929) < eps,
			"Incorrect K at 0, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 4) -
				0.9787589339778296) < eps,
			"Incorrect K at 0, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 5) -
				0.9904431525124753) < eps,
			"Incorrect K at 0, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 6) -
				0.9812018448276493) < eps,
			"Incorrect K at 0, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 7) -
				0.9846502712939124) < eps,
			"Incorrect K at 0, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 8) -
				0.9852360775152181) < eps,
			"Incorrect K at 0, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 0, 9) -
				0.9726553489464377) < eps,
			"Incorrect K at 0, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 0) -
				0.9946686121177555) < eps,
			"Incorrect K at 1, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 1) -
				0.9972853441139281) < eps,
			"Incorrect K at 1, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 2) -
				0.9737099505662699) < eps,
			"Incorrect K at 1, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 3) -
				0.9873570938969980) < eps,
			"Incorrect K at 1, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 4) -
				0.9872109867697908) < eps,
			"Incorrect K at 1, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 5) -
				0.9849733935382557) < eps,
			"Incorrect K at 1, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 6) -
				0.9885619953841569) < eps,
			"Incorrect K at 1, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 7) -
				0.9736717805509784) < eps,
			"Incorrect K at 1, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 8) -
				0.9841030912139471) < eps,
			"Incorrect K at 1, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 1, 9) -
				0.9651822304880787) < eps,
			"Incorrect K at 1, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 0) -
				0.9578937155783831) < eps,
			"Incorrect K at 2, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 1) -
				0.9737099505662699) < eps,
			"Incorrect K at 2, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 2) -
				0.9574301966133943) < eps,
			"Incorrect K at 2, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 3) -
				0.9611091259460700) < eps,
			"Incorrect K at 2, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 4) -
				0.9709023602232207) < eps,
			"Incorrect K at 2, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 5) -
				0.9315907625390656) < eps,
			"Incorrect K at 2, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 6) -
				0.9705232203974647) < eps,
			"Incorrect K at 2, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 7) -
				0.9348012913572012) < eps,
			"Incorrect K at 2, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 8) -
				0.9594537644613513) < eps,
			"Incorrect K at 2, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 2, 9) -
				0.9323746182879774) < eps,
			"Incorrect K at 2, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 0) -
				0.9892315762620929) < eps,
			"Incorrect K at 3, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 1) -
				0.9873570938969980) < eps,
			"Incorrect K at 3, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 2) -
				0.9611091259460700) < eps,
			"Incorrect K at 3, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 3) -
				0.9951110594420507) < eps,
			"Incorrect K at 3, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 4) -
				0.9877282227805176) < eps,
			"Incorrect K at 3, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 5) -
				0.9865497273302540) < eps,
			"Incorrect K at 3, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 6) -
				0.9876865775961812) < eps,
			"Incorrect K at 3, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 7) -
				0.9906424333748519) < eps,
			"Incorrect K at 3, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 8) -
				0.9928175000327435) < eps,
			"Incorrect K at 3, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 3, 9) -
				0.9758002897316713) < eps,
			"Incorrect K at 3, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 0) -
				0.9787589339778296) < eps,
			"Incorrect K at 4, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 1) -
				0.9872109867697908) < eps,
			"Incorrect K at 4, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 2) -
				0.9709023602232207) < eps,
			"Incorrect K at 4, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 3) -
				0.9877282227805176) < eps,
			"Incorrect K at 4, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 4) -
				0.9882189022474748) < eps,
			"Incorrect K at 4, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 5) -
				0.9599174207778631) < eps,
			"Incorrect K at 4, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 6) -
				0.9877988987284404) < eps,
			"Incorrect K at 4, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 7) -
				0.9700057126716283) < eps,
			"Incorrect K at 4, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 8) -
				0.9852740882958205) < eps,
			"Incorrect K at 4, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 4, 9) -
				0.9555553363633501) < eps,
			"Incorrect K at 4, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 0) -
				0.9904431525124753) < eps,
			"Incorrect K at 5, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 1) -
				0.9849733935382557) < eps,
			"Incorrect K at 5, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 2) -
				0.9315907625390656) < eps,
			"Incorrect K at 5, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 3) -
				0.9865497273302540) < eps,
			"Incorrect K at 5, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 4) -
				0.9599174207778631) < eps,
			"Incorrect K at 5, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 5) -
				0.9902416959670991) < eps,
			"Incorrect K at 5, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 6) -
				0.9639775155470058) < eps,
			"Incorrect K at 5, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 7) -
				0.9867152808385368) < eps,
			"Incorrect K at 5, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 8) -
				0.9813171680326727) < eps,
			"Incorrect K at 5, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 5, 9) -
				0.9727966479866800) < eps,
			"Incorrect K at 5, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 0) -
				0.9812018448276493) < eps,
			"Incorrect K at 6, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 1) -
				0.9885619953841569) < eps,
			"Incorrect K at 6, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 2) -
				0.9705232203974647) < eps,
			"Incorrect K at 6, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 3) -
				0.9876865775961812) < eps,
			"Incorrect K at 6, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 4) -
				0.9877988987284404) < eps,
			"Incorrect K at 6, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 5) -
				0.9639775155470058) < eps,
			"Incorrect K at 6, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 6) -
				0.9875439243202984) < eps,
			"Incorrect K at 6, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 7) -
				0.9708912767839261) < eps,
			"Incorrect K at 6, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 8) -
				0.9851182379707764) < eps,
			"Incorrect K at 6, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 6, 9) -
				0.9568814495759891) < eps,
			"Incorrect K at 6, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 0) -
				0.9846502712939124) < eps,
			"Incorrect K at 7, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 1) -
				0.9736717805509784) < eps,
			"Incorrect K at 7, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 2) -
				0.9348012913572012) < eps,
			"Incorrect K at 7, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 3) -
				0.9906424333748519) < eps,
			"Incorrect K at 7, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 4) -
				0.9700057126716283) < eps,
			"Incorrect K at 7, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 5) -
				0.9867152808385368) < eps,
			"Incorrect K at 7, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 6) -
				0.9708912767839261) < eps,
			"Incorrect K at 7, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 7) -
				0.9889671313757883) < eps,
			"Incorrect K at 7, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 8) -
				0.9867135908997623) < eps,
			"Incorrect K at 7, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 7, 9) -
				0.9731498666169363) < eps,
			"Incorrect K at 7, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 0) -
				0.9852360775152181) < eps,
			"Incorrect K at 8, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 1) -
				0.9841030912139471) < eps,
			"Incorrect K at 8, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 2) -
				0.9594537644613513) < eps,
			"Incorrect K at 8, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 3) -
				0.9928175000327435) < eps,
			"Incorrect K at 8, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 4) -
				0.9852740882958205) < eps,
			"Incorrect K at 8, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 5) -
				0.9813171680326727) < eps,
			"Incorrect K at 8, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 6) -
				0.9851182379707764) < eps,
			"Incorrect K at 8, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 7) -
				0.9867135908997623) < eps,
			"Incorrect K at 8, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 8) -
				0.9900919372870961) < eps,
			"Incorrect K at 8, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 8, 9) -
				0.9706093824836438) < eps,
			"Incorrect K at 8, 9");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 0) -
				0.9726553489464377) < eps,
			"Incorrect K at 9, 0");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 1) -
				0.9651822304880787) < eps,
			"Incorrect K at 9, 1");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 2) -
				0.9323746182879774) < eps,
			"Incorrect K at 9, 2");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 3) -
				0.9758002897316713) < eps,
			"Incorrect K at 9, 3");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 4) -
				0.9555553363633501) < eps,
			"Incorrect K at 9, 4");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 5) -
				0.9727966479866800) < eps,
			"Incorrect K at 9, 5");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 6) -
				0.9568814495759891) < eps,
			"Incorrect K at 9, 6");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 7) -
				0.9731498666169363) < eps,
			"Incorrect K at 9, 7");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 8) -
				0.9706093824836438) < eps,
			"Incorrect K at 9, 8");
	mu_assert(fabs(matrix_get(K, data->n, data->n, 9, 9) -
				0.9572935857028608) < eps,
			"Incorrect K at 9, 9");

	// end test code //

	free(K);
	gensvm_free_data(data);
	gensvm_free_model(model);

	return NULL;
}

char *test_kernel_eigendecomp()
{
	int n = 10;
	double *K = Calloc(double, n*n);

	matrix_set(K, n, n, 0, 0, 1.8957230074596658);
	matrix_set(K, n, n, 0, 1, 1.5788353939003679);
	matrix_set(K, n, n, 0, 2, 1.3341924479973508);
	matrix_set(K, n, n, 0, 3, 1.9135891286401130);
	matrix_set(K, n, n, 0, 4, 1.0234654124933122);
	matrix_set(K, n, n, 0, 5, 1.2788676188422881);
	matrix_set(K, n, n, 0, 6, 1.5231707888435004);
	matrix_set(K, n, n, 0, 7, 1.7650988656771305);
	matrix_set(K, n, n, 0, 8, 1.6565264972793248);
	matrix_set(K, n, n, 0, 9, 1.7720466624967910);
	matrix_set(K, n, n, 1, 0, 1.5788353939003679);
	matrix_set(K, n, n, 1, 1, 3.3160214905692085);
	matrix_set(K, n, n, 1, 2, 1.9970415315746892);
	matrix_set(K, n, n, 1, 3, 2.2818422166445553);
	matrix_set(K, n, n, 1, 4, 2.2552591311973780);
	matrix_set(K, n, n, 1, 5, 2.6229134325776573);
	matrix_set(K, n, n, 1, 6, 2.1715097610253036);
	matrix_set(K, n, n, 1, 7, 1.6910997290719514);
	matrix_set(K, n, n, 1, 8, 2.6386724572773361);
	matrix_set(K, n, n, 1, 9, 2.6587015780393881);
	matrix_set(K, n, n, 2, 0, 1.3341924479973508);
	matrix_set(K, n, n, 2, 1, 1.9970415315746892);
	matrix_set(K, n, n, 2, 2, 2.9794552376145154);
	matrix_set(K, n, n, 2, 3, 2.6421606833813742);
	matrix_set(K, n, n, 2, 4, 2.2526586210360087);
	matrix_set(K, n, n, 2, 5, 2.0022301831640106);
	matrix_set(K, n, n, 2, 6, 2.5890262336485397);
	matrix_set(K, n, n, 2, 7, 2.2391395202916926);
	matrix_set(K, n, n, 2, 8, 2.4893629764352507);
	matrix_set(K, n, n, 2, 9, 2.7339053828614057);
	matrix_set(K, n, n, 3, 0, 1.9135891286401130);
	matrix_set(K, n, n, 3, 1, 2.2818422166445553);
	matrix_set(K, n, n, 3, 2, 2.6421606833813742);
	matrix_set(K, n, n, 3, 3, 3.4991899582151325);
	matrix_set(K, n, n, 3, 4, 1.8812151717558272);
	matrix_set(K, n, n, 3, 5, 2.6296896695972061);
	matrix_set(K, n, n, 3, 6, 2.8161237364134410);
	matrix_set(K, n, n, 3, 7, 2.7729181198442641);
	matrix_set(K, n, n, 3, 8, 2.8241735649875288);
	matrix_set(K, n, n, 3, 9, 2.9235411206761861);
	matrix_set(K, n, n, 4, 0, 1.0234654124933122);
	matrix_set(K, n, n, 4, 1, 2.2552591311973780);
	matrix_set(K, n, n, 4, 2, 2.2526586210360087);
	matrix_set(K, n, n, 4, 3, 1.8812151717558272);
	matrix_set(K, n, n, 4, 4, 2.5093034717737615);
	matrix_set(K, n, n, 4, 5, 1.8901064622829966);
	matrix_set(K, n, n, 4, 6, 1.8736444909319021);
	matrix_set(K, n, n, 4, 7, 1.4654421216986091);
	matrix_set(K, n, n, 4, 8, 2.2618348136979898);
	matrix_set(K, n, n, 4, 9, 1.9853102751114591);
	matrix_set(K, n, n, 5, 0, 1.2788676188422881);
	matrix_set(K, n, n, 5, 1, 2.6229134325776573);
	matrix_set(K, n, n, 5, 2, 2.0022301831640106);
	matrix_set(K, n, n, 5, 3, 2.6296896695972061);
	matrix_set(K, n, n, 5, 4, 1.8901064622829966);
	matrix_set(K, n, n, 5, 5, 2.8486874970858160);
	matrix_set(K, n, n, 5, 6, 2.1833913326358250);
	matrix_set(K, n, n, 5, 7, 1.6103467044130788);
	matrix_set(K, n, n, 5, 8, 2.5103824469216467);
	matrix_set(K, n, n, 5, 9, 2.5894371939176293);
	matrix_set(K, n, n, 6, 0, 1.5231707888435004);
	matrix_set(K, n, n, 6, 1, 2.1715097610253036);
	matrix_set(K, n, n, 6, 2, 2.5890262336485397);
	matrix_set(K, n, n, 6, 3, 2.8161237364134410);
	matrix_set(K, n, n, 6, 4, 1.8736444909319021);
	matrix_set(K, n, n, 6, 5, 2.1833913326358250);
	matrix_set(K, n, n, 6, 6, 2.9754371323085365);
	matrix_set(K, n, n, 6, 7, 2.0569948157192495);
	matrix_set(K, n, n, 6, 8, 2.7702792216087251);
	matrix_set(K, n, n, 6, 9, 2.6955896898884357);
	matrix_set(K, n, n, 7, 0, 1.7650988656771305);
	matrix_set(K, n, n, 7, 1, 1.6910997290719514);
	matrix_set(K, n, n, 7, 2, 2.2391395202916926);
	matrix_set(K, n, n, 7, 3, 2.7729181198442641);
	matrix_set(K, n, n, 7, 4, 1.4654421216986091);
	matrix_set(K, n, n, 7, 5, 1.6103467044130788);
	matrix_set(K, n, n, 7, 6, 2.0569948157192495);
	matrix_set(K, n, n, 7, 7, 2.8491469863450520);
	matrix_set(K, n, n, 7, 8, 1.9828570918824218);
	matrix_set(K, n, n, 7, 9, 2.3284525338942448);
	matrix_set(K, n, n, 8, 0, 1.6565264972793248);
	matrix_set(K, n, n, 8, 1, 2.6386724572773361);
	matrix_set(K, n, n, 8, 2, 2.4893629764352507);
	matrix_set(K, n, n, 8, 3, 2.8241735649875288);
	matrix_set(K, n, n, 8, 4, 2.2618348136979898);
	matrix_set(K, n, n, 8, 5, 2.5103824469216467);
	matrix_set(K, n, n, 8, 6, 2.7702792216087251);
	matrix_set(K, n, n, 8, 7, 1.9828570918824218);
	matrix_set(K, n, n, 8, 8, 3.1036474712902948);
	matrix_set(K, n, n, 8, 9, 2.7267087162991572);
	matrix_set(K, n, n, 9, 0, 1.7720466624967910);
	matrix_set(K, n, n, 9, 1, 2.6587015780393881);
	matrix_set(K, n, n, 9, 2, 2.7339053828614057);
	matrix_set(K, n, n, 9, 3, 2.9235411206761861);
	matrix_set(K, n, n, 9, 4, 1.9853102751114591);
	matrix_set(K, n, n, 9, 5, 2.5894371939176293);
	matrix_set(K, n, n, 9, 6, 2.6955896898884357);
	matrix_set(K, n, n, 9, 7, 2.3284525338942448);
	matrix_set(K, n, n, 9, 8, 2.7267087162991572);
	matrix_set(K, n, n, 9, 9, 3.6481211298649119);

	// start test code //
	double *P = NULL;
	double *Sigma = NULL;
	long r = gensvm_kernel_eigendecomp(K, n, 1e-2, &P, &Sigma);
	double eps = 1e-13;

	mu_assert(r == 7, "Incorrect number of eigenvalues kept");

	// Note: to overcome sign variability in the eigenvectors, we take the
	// absolute value of the elements of P and the expected outcome.

	mu_assert(fabs(fabs(matrix_get(P, n, r, 0, 0)) -
				fabs(0.2153575812884510)) < eps,
			"Incorrect P at 0, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 0, 1)) -
				fabs(0.2745973925101974)) < eps,
			"Incorrect P at 0, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 0, 2)) -
				fabs(0.4550895791914591)) < eps,
			"Incorrect P at 0, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 0, 3)) -
				fabs(-0.3884421314028316)) < eps,
			"Incorrect P at 0, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 0, 4)) -
				fabs(0.1160315888082866)) < eps,
			"Incorrect P at 0, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 0, 5)) -
				fabs(0.4302035203307312)) < eps,
			"Incorrect P at 0, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 0, 6)) -
				fabs(-0.4669338527728559)) < eps,
			"Incorrect P at 0, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 1, 0)) -
				fabs(0.3217489710376946)) < eps,
			"Incorrect P at 1, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 1, 1)) -
				fabs(-0.4945822605094066)) < eps,
			"Incorrect P at 1, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 1, 2)) -
				fabs(0.3769753694343475)) < eps,
			"Incorrect P at 1, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 1, 3)) -
				fabs(-0.2972948635736606)) < eps,
			"Incorrect P at 1, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 1, 4)) -
				fabs(-0.1011097319500622)) < eps,
			"Incorrect P at 1, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 1, 5)) -
				fabs(0.0306004375001570)) < eps,
			"Incorrect P at 1, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 1, 6)) -
				fabs(0.5322704452300250)) < eps,
			"Incorrect P at 1, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 2, 0)) -
				fabs(0.3239269538650612)) < eps,
			"Incorrect P at 2, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 2, 1)) -
				fabs(0.1131046585105233)) < eps,
			"Incorrect P at 2, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 2, 2)) -
				fabs(-0.5683419560926302)) < eps,
			"Incorrect P at 2, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 2, 3)) -
				fabs(-0.0566442945574068)) < eps,
			"Incorrect P at 2, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 2, 4)) -
				fabs(-0.2195800248661894)) < eps,
			"Incorrect P at 2, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 2, 5)) -
				fabs(-0.0120101816760300)) < eps,
			"Incorrect P at 2, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 2, 6)) -
				fabs(-0.0868590923278641)) < eps,
			"Incorrect P at 2, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 3, 0)) -
				fabs(0.3643049297817086)) < eps,
			"Incorrect P at 3, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 3, 1)) -
				fabs(0.3137060295518203)) < eps,
			"Incorrect P at 3, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 3, 2)) -
				fabs(0.1059271459195730)) < eps,
			"Incorrect P at 3, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 3, 3)) -
				fabs(0.2260527473443755)) < eps,
			"Incorrect P at 3, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 3, 4)) -
				fabs(0.2873469240341838)) < eps,
			"Incorrect P at 3, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 3, 5)) -
				fabs(-0.3577074997001115)) < eps,
			"Incorrect P at 3, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 3, 6)) -
				fabs(-0.1176008886950781)) < eps,
			"Incorrect P at 3, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 4, 0)) -
				fabs(0.2685490754867729)) < eps,
			"Incorrect P at 4, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 4, 1)) -
				fabs(-0.3602884322676240)) < eps,
			"Incorrect P at 4, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 4, 2)) -
				fabs(-0.4208535241232802)) < eps,
			"Incorrect P at 4, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 4, 3)) -
				fabs(-0.5004289581783992)) < eps,
			"Incorrect P at 4, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 4, 4)) -
				fabs(0.0581201143184524)) < eps,
			"Incorrect P at 4, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 4, 5)) -
				fabs(-0.1560533086067875)) < eps,
			"Incorrect P at 4, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 4, 6)) -
				fabs(-0.2808692440831040)) < eps,
			"Incorrect P at 4, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 5, 0)) -
				fabs(0.3097165272244328)) < eps,
			"Incorrect P at 5, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 5, 1)) -
				fabs(-0.3150295049924597)) < eps,
			"Incorrect P at 5, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 5, 2)) -
				fabs(0.2661795640562757)) < eps,
			"Incorrect P at 5, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 5, 3)) -
				fabs(0.3687999410293836)) < eps,
			"Incorrect P at 5, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 5, 4)) -
				fabs(0.0756501054032381)) < eps,
			"Incorrect P at 5, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 5, 5)) -
				fabs(-0.4702905504198338)) < eps,
			"Incorrect P at 5, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 5, 6)) -
				fabs(-0.2502907081561355)) < eps,
			"Incorrect P at 5, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 6, 0)) -
				fabs(0.3303343184191819)) < eps,
			"Incorrect P at 6, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 6, 1)) -
				fabs(0.0877791634967746)) < eps,
			"Incorrect P at 6, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 6, 2)) -
				fabs(-0.2291903024112904)) < eps,
			"Incorrect P at 6, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 6, 3)) -
				fabs(0.3378230216485225)) < eps,
			"Incorrect P at 6, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 6, 4)) -
				fabs(0.2813505143218815)) < eps,
			"Incorrect P at 6, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 6, 5)) -
				fabs(0.4390481386320075)) < eps,
			"Incorrect P at 6, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 6, 6)) -
				fabs(0.3960648555497706)) < eps,
			"Incorrect P at 6, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 7, 0)) -
				fabs(0.2867008390424503)) < eps,
			"Incorrect P at 7, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 7, 1)) -
				fabs(0.5529782699385689)) < eps,
			"Incorrect P at 7, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 7, 2)) -
				fabs(0.0595354369281301)) < eps,
			"Incorrect P at 7, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 7, 3)) -
				fabs(-0.3535615166878320)) < eps,
			"Incorrect P at 7, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 7, 4)) -
				fabs(-0.0706307412848860)) < eps,
			"Incorrect P at 7, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 7, 5)) -
				fabs(-0.3515408303354066)) < eps,
			"Incorrect P at 7, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 7, 6)) -
				fabs(0.3840483694579672)) < eps,
			"Incorrect P at 7, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 8, 0)) -
				fabs(0.3474605686018905)) < eps,
			"Incorrect P at 8, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 8, 1)) -
				fabs(-0.1590687116733123)) < eps,
			"Incorrect P at 8, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 8, 2)) -
				fabs(-0.0565892834229501)) < eps,
			"Incorrect P at 8, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 8, 3)) -
				fabs(0.0725995550169575)) < eps,
			"Incorrect P at 8, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 8, 4)) -
				fabs(0.4060731268950495)) < eps,
			"Incorrect P at 8, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 8, 5)) -
				fabs(0.2797250817800897)) < eps,
			"Incorrect P at 8, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 8, 6)) -
				fabs(-0.0976363581213828)) < eps,
			"Incorrect P at 8, 6");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 9, 0)) -
				fabs(0.3638159501275298)) < eps,
			"Incorrect P at 9, 0");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 9, 1)) -
				fabs(0.0305969903511121)) < eps,
			"Incorrect P at 9, 1");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 9, 2)) -
				fabs(0.0964657835567429)) < eps,
			"Incorrect P at 9, 2");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 9, 3)) -
				fabs(0.2749121311297592)) < eps,
			"Incorrect P at 9, 3");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 9, 4)) -
				fabs(-0.7664116168127906)) < eps,
			"Incorrect P at 9, 4");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 9, 5)) -
				fabs(0.2139473364333364)) < eps,
			"Incorrect P at 9, 5");
	mu_assert(fabs(fabs(matrix_get(P, n, r, 9, 6)) -
				fabs(-0.1478517167777606)) < eps,
			"Incorrect P at 9, 6");

	eps = 1e-13;

	mu_assert(fabs(Sigma[0] - 4.7923242920575353) < eps,
			"Incorrect Sigma at 0");
	mu_assert(fabs(Sigma[1] - 1.5023267732150303) < eps,
			"Incorrect Sigma at 1");
	mu_assert(fabs(Sigma[2] - 1.1906890711161726) < eps,
			"Incorrect Sigma at 2");
	mu_assert(fabs(Sigma[3] - 1.0037677343120393) < eps,
			"Incorrect Sigma at 3");
	mu_assert(fabs(Sigma[4] - 0.8899042921295091) < eps,
			"Incorrect Sigma at 4");
	mu_assert(fabs(Sigma[5] - 0.8251500911319253) < eps,
			"Incorrect Sigma at 5");
	mu_assert(fabs(Sigma[6] - 0.5394180349552358) < eps,
			"Incorrect Sigma at 6");

	// end test code //

	free(K);
	free(P);
	free(Sigma);

	return NULL;
}

char *test_kernel_cross_rbf()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data_1 = gensvm_init_data();
	struct GenData *data_2 = gensvm_init_data();

	// setup //
	data_1->n = 10;
	data_1->m = 3;
	data_1->RAW = Calloc(double, data_1->n * (data_1->m + 1));

	data_2->n = 5;
	data_2->m = 3;
	data_2->RAW = Calloc(double, data_2->n * (data_2->m + 1));

	model->n = 10;
	model->m = 3;
	model->kerneltype = K_RBF;
	model->gamma = 0.348;

	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 3, 0.0884616753393881);

	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 1, 0.8233234072519983);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 2, 0.3267543833513200);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 3, 0.2728942849022845);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 1, 0.7956168453294307);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 2, 0.3693175185473680);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 3, 0.2665038412777220);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 1, 0.5777227081256917);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 2, 0.2456426390463990);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 3, 0.1126634767453339);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 1, 0.4426030703804438);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 2, 0.1102697774064020);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 3, 0.9577384085295405);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 1, 0.7701104553132680);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 2, 0.8659836346403005);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 3, 0.4421972387164297);

	// start test code //
	double eps = 1e-14;
	double *K2 = gensvm_kernel_cross(model, data_1, data_2);
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 0) -
				0.9807518230471629) < eps,
			"Incorrect K2 at 0, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 1) -
				0.8852741170384059) < eps,
			"Incorrect K2 at 0, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 2) -
				0.7635716665765656) < eps,
			"Incorrect K2 at 0, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 3) -
				0.7492080432261149) < eps,
			"Incorrect K2 at 0, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 4) -
				0.7170903046043821) < eps,
			"Incorrect K2 at 0, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 5) -
				0.9391346724948219) < eps,
			"Incorrect K2 at 0, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 6) -
				0.7511075325872592) < eps,
			"Incorrect K2 at 0, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 7) -
				0.8068193111978309) < eps,
			"Incorrect K2 at 0, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 8) -
				0.7712432709418185) < eps,
			"Incorrect K2 at 0, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 9) -
				0.8700312205062914) < eps,
			"Incorrect K2 at 0, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 0) -
				0.9841138078975342) < eps,
			"Incorrect K2 at 1, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 1) -
				0.8790397003296203) < eps,
			"Incorrect K2 at 1, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 2) -
				0.7669459160670393) < eps,
			"Incorrect K2 at 1, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 3) -
				0.7703671633763146) < eps,
			"Incorrect K2 at 1, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 4) -
				0.7272525607549757) < eps,
			"Incorrect K2 at 1, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 5) -
				0.9511025236778456) < eps,
			"Incorrect K2 at 1, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 6) -
				0.7605364853756941) < eps,
			"Incorrect K2 at 1, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 7) -
				0.8282285706223327) < eps,
			"Incorrect K2 at 1, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 8) -
				0.7905551836844474) < eps,
			"Incorrect K2 at 1, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 9) -
				0.8829411638604122) < eps,
			"Incorrect K2 at 1, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 0) -
				0.9259609685011719) < eps,
			"Incorrect K2 at 2, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 1) -
				0.8145263140370310) < eps,
			"Incorrect K2 at 2, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 2) -
				0.8395046061962322) < eps,
			"Incorrect K2 at 2, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 3) -
				0.7622039915518654) < eps,
			"Incorrect K2 at 2, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 4) -
				0.7511426474507734) < eps,
			"Incorrect K2 at 2, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 5) -
				0.9415384710796756) < eps,
			"Incorrect K2 at 2, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 6) -
				0.7783683734706774) < eps,
			"Incorrect K2 at 2, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 7) -
				0.8503122106685652) < eps,
			"Incorrect K2 at 2, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 8) -
				0.8011956519011766) < eps,
			"Incorrect K2 at 2, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 9) -
				0.9458326999989332) < eps,
			"Incorrect K2 at 2, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 0) -
				0.8296567517703916) < eps,
			"Incorrect K2 at 3, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 1) -
				0.9519346106202000) < eps,
			"Incorrect K2 at 3, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 2) -
				0.8748910369012700) < eps,
			"Incorrect K2 at 3, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 3) -
				0.7304409952010390) < eps,
			"Incorrect K2 at 3, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 4) -
				0.9040565707355339) < eps,
			"Incorrect K2 at 3, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 5) -
				0.6649291594794616) < eps,
			"Incorrect K2 at 3, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 6) -
				0.9229382678947169) < eps,
			"Incorrect K2 at 3, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 7) -
				0.6306350314426549) < eps,
			"Incorrect K2 at 3, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 8) -
				0.7718806233578666) < eps,
			"Incorrect K2 at 3, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 9) -
				0.7285359056939609) < eps,
			"Incorrect K2 at 3, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 0) -
				0.9509298370980613) < eps,
			"Incorrect K2 at 4, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 1) -
				0.8040728104591074) < eps,
			"Incorrect K2 at 4, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 2) -
				0.6448137876591025) < eps,
			"Incorrect K2 at 4, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 3) -
				0.8810145544703872) < eps,
			"Incorrect K2 at 4, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 4) -
				0.7167145071888406) < eps,
			"Incorrect K2 at 4, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 5) -
				0.9268433969926227) < eps,
			"Incorrect K2 at 4, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 6) -
				0.7452783549758193) < eps,
			"Incorrect K2 at 4, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 7) -
				0.8788513369285231) < eps,
			"Incorrect K2 at 4, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 8) -
				0.8577568482722602) < eps,
			"Incorrect K2 at 4, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 9) -
				0.8117847908337736) < eps,
			"Incorrect K2 at 4, 9");
	free(K2);

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data_1);
	gensvm_free_data(data_2);

	return NULL;
}

char *test_kernel_cross_poly()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data_1 = gensvm_init_data();
	struct GenData *data_2 = gensvm_init_data();

	// setup //
	data_1->n = 10;
	data_1->m = 3;
	data_1->RAW = Calloc(double, data_1->n * (data_1->m + 1));

	data_2->n = 5;
	data_2->m = 3;
	data_2->RAW = Calloc(double, data_2->n * (data_2->m + 1));

	model->n = 10;
	model->m = 3;
	model->kerneltype = K_POLY;
	model->gamma = 1.5;
	model->coef = 3.0;
	model->degree = 1.78;

	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 3, 0.0884616753393881);

	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 1, 0.8233234072519983);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 2, 0.3267543833513200);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 3, 0.2728942849022845);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 1, 0.7956168453294307);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 2, 0.3693175185473680);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 3, 0.2665038412777220);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 1, 0.5777227081256917);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 2, 0.2456426390463990);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 3, 0.1126634767453339);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 1, 0.4426030703804438);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 2, 0.1102697774064020);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 3, 0.9577384085295405);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 1, 0.7701104553132680);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 2, 0.8659836346403005);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 3, 0.4421972387164297);

	// start test code //
	double eps = 1e-14;
	double *K2 = gensvm_kernel_cross(model, data_1, data_2);
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 0) -
				14.0660499715713421) < eps,
			"Incorrect K2 at 0, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 1) -
				14.0798736090876808) < eps,
			"Incorrect K2 at 0, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 2) -
				8.1700400917705256) < eps,
			"Incorrect K2 at 0, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 3) -
				11.1676455855990469) < eps,
			"Incorrect K2 at 0, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 4) -
				9.3728783610105122) < eps,
			"Incorrect K2 at 0, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 5) -
				12.6182273329088108) < eps,
			"Incorrect K2 at 0, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 6) -
				9.7694076615604182) < eps,
			"Incorrect K2 at 0, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 7) -
				10.7135530249071351) < eps,
			"Incorrect K2 at 0, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 8) -
				10.3924703115971937) < eps,
			"Incorrect K2 at 0, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 9) -
				9.4643617813593419) < eps,
			"Incorrect K2 at 0, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 0) -
				14.0284757046542925) < eps,
			"Incorrect K2 at 1, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 1) -
				13.9144650311002316) < eps,
			"Incorrect K2 at 1, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 2) -
				8.1499265608165583) < eps,
			"Incorrect K2 at 1, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 3) -
				11.4045391037528390) < eps,
			"Incorrect K2 at 1, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 4) -
				9.4504735861367735) < eps,
			"Incorrect K2 at 1, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 5) -
				12.6901332253769485) < eps,
			"Incorrect K2 at 1, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 6) -
				9.8317381138302018) < eps,
			"Incorrect K2 at 1, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 7) -
				10.9280234697682062) < eps,
			"Incorrect K2 at 1, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 8) -
				10.5883823060327860) < eps,
			"Incorrect K2 at 1, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 9) -
				9.5490652847135138) < eps,
			"Incorrect K2 at 1, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 0) -
				11.5485754684510606) < eps,
			"Incorrect K2 at 2, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 1) -
				11.2735171075936282) < eps,
			"Incorrect K2 at 2, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 2) -
				7.5769874635472778) < eps,
			"Incorrect K2 at 2, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 3) -
				9.6665605291258849) < eps,
			"Incorrect K2 at 2, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 4) -
				8.2678047879187648) < eps,
			"Incorrect K2 at 2, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 5) -
				10.8709836661754018) < eps,
			"Incorrect K2 at 2, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 6) -
				8.5376169230827301) < eps,
			"Incorrect K2 at 2, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 7) -
				9.6025396489681700) < eps,
			"Incorrect K2 at 2, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 8) -
				9.1503562342971971) < eps,
			"Incorrect K2 at 2, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 9) -
				8.7190361363030764) < eps,
			"Incorrect K2 at 2, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 0) -
				13.1667427919927871) < eps,
			"Incorrect K2 at 3, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 1) -
				16.1629665429693823) < eps,
			"Incorrect K2 at 3, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 2) -
				10.5015908640040099) < eps,
			"Incorrect K2 at 3, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 3) -
				11.9214318530013443) < eps,
			"Incorrect K2 at 3, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 4) -
				12.9437638270688211) < eps,
			"Incorrect K2 at 3, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 5) -
				9.8267375915304900) < eps,
			"Incorrect K2 at 3, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 6) -
				13.0942580006333298) < eps,
			"Incorrect K2 at 3, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 7) -
				9.1036164336081900) < eps,
			"Incorrect K2 at 3, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 8) -
				11.4137186988419348) < eps,
			"Incorrect K2 at 3, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 9) -
				8.6105081033670565) < eps,
			"Incorrect K2 at 3, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 0) -
				16.6753988505663493) < eps,
			"Incorrect K2 at 4, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 1) -
				15.8396512236024520) < eps,
			"Incorrect K2 at 4, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 2) -
				8.8303654890651160) < eps,
			"Incorrect K2 at 4, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 3) -
				15.9428461498604186) < eps,
			"Incorrect K2 at 4, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 4) -
				11.9205101781336413) < eps,
			"Incorrect K2 at 4, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 5) -
				15.3360611842713368) < eps,
			"Incorrect K2 at 4, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 6) -
				12.2768145090113094) < eps,
			"Incorrect K2 at 4, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 7) -
				14.4512803240791001) < eps,
			"Incorrect K2 at 4, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 8) -
				14.3401441810752441) < eps,
			"Incorrect K2 at 4, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 9) -
				11.2489550339576638) < eps,
			"Incorrect K2 at 4, 9");
	free(K2);

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data_1);
	gensvm_free_data(data_2);

	return NULL;
}

char *test_kernel_cross_sigmoid()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data_1 = gensvm_init_data();
	struct GenData *data_2 = gensvm_init_data();

	// setup //
	data_1->n = 10;
	data_1->m = 3;
	data_1->RAW = Calloc(double, data_1->n * (data_1->m + 1));

	data_2->n = 5;
	data_2->m = 3;
	data_2->RAW = Calloc(double, data_2->n * (data_2->m + 1));

	model->n = 10;
	model->m = 3;
	model->kerneltype = K_SIGMOID;
	model->gamma = 1.23;
	model->coef = 1.6;

	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 0, 1.0000000000000000);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data_1->RAW, data_1->n, data_1->m+1, 9, 3, 0.0884616753393881);

	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 1, 0.8233234072519983);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 2, 0.3267543833513200);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 0, 3, 0.2728942849022845);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 1, 0.7956168453294307);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 2, 0.3693175185473680);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 1, 3, 0.2665038412777220);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 1, 0.5777227081256917);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 2, 0.2456426390463990);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 2, 3, 0.1126634767453339);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 1, 0.4426030703804438);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 2, 0.1102697774064020);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 3, 3, 0.9577384085295405);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 1, 0.7701104553132680);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 2, 0.8659836346403005);
	matrix_set(data_2->RAW, data_2->n, data_2->m+1, 4, 3, 0.4421972387164297);

	// start test code //
	double eps = 1e-14;
	double *K2 = gensvm_kernel_cross(model, data_1, data_2);

	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 0) -
				0.9920395702373311) < eps,
			"Incorrect K2 at 0, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 1) -
				0.9920712050042086) < eps,
			"Incorrect K2 at 0, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 2) -
				0.9476993343669046) < eps,
			"Incorrect K2 at 0, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 3) -
				0.9809051895947648) < eps,
			"Incorrect K2 at 0, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 4) -
				0.9656014407964381) < eps,
			"Incorrect K2 at 0, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 5) -
				0.9878041245323028) < eps,
			"Incorrect K2 at 0, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 6) -
				0.9699033492645315) < eps,
			"Incorrect K2 at 0, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 7) -
				0.9779225156794014) < eps,
			"Incorrect K2 at 0, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 8) -
				0.9754999374304111) < eps,
			"Incorrect K2 at 0, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 0, 9) -
				0.9666518758877651) < eps,
			"Incorrect K2 at 0, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 0) -
				0.9919528767612884) < eps,
			"Incorrect K2 at 1, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 1) -
				0.9916833786531655) < eps,
			"Incorrect K2 at 1, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 2) -
				0.9473218250023265) < eps,
			"Incorrect K2 at 1, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 3) -
				0.9822808818850081) < eps,
			"Incorrect K2 at 1, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 4) -
				0.9664947467723793) < eps,
			"Incorrect K2 at 1, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 5) -
				0.9880653879683305) < eps,
			"Incorrect K2 at 1, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 6) -
				0.9705231426655428) < eps,
			"Incorrect K2 at 1, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 7) -
				0.9793914678629526) < eps,
			"Incorrect K2 at 1, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 8) -
				0.9770114192540243) < eps,
			"Incorrect K2 at 1, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 1, 9) -
				0.9675925213675517) < eps,
			"Incorrect K2 at 1, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 0) -
				0.9830631959900493) < eps,
			"Incorrect K2 at 2, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 1) -
				0.9815341888355711) < eps,
			"Incorrect K2 at 2, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 2) -
				0.9351678264557334) < eps,
			"Incorrect K2 at 2, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 3) -
				0.9688485293248640) < eps,
			"Incorrect K2 at 2, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 4) -
				0.9494915246914850) < eps,
			"Incorrect K2 at 2, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 5) -
				0.9790117279097025) < eps,
			"Incorrect K2 at 2, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 6) -
				0.9540901254104273) < eps,
			"Incorrect K2 at 2, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 7) -
				0.9681710508625928) < eps,
			"Incorrect K2 at 2, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 8) -
				0.9628881304379768) < eps,
			"Incorrect K2 at 2, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 2, 9) -
				0.9569175771195083) < eps,
			"Incorrect K2 at 2, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 0) -
				0.9896480192770086) < eps,
			"Incorrect K2 at 3, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 1) -
				0.9955708509948572) < eps,
			"Incorrect K2 at 3, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 2) -
				0.9763550511314351) < eps,
			"Incorrect K2 at 3, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 3) -
				0.9849157087004718) < eps,
			"Incorrect K2 at 3, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 4) -
				0.9889385756779369) < eps,
			"Incorrect K2 at 3, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 5) -
				0.9704739495606978) < eps,
			"Incorrect K2 at 3, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 6) -
				0.9894230831842615) < eps,
			"Incorrect K2 at 3, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 7) -
				0.9622882117865678) < eps,
			"Incorrect K2 at 3, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 8) -
				0.9823319222305995) < eps,
			"Incorrect K2 at 3, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 3, 9) -
				0.9552504472873038) < eps,
			"Incorrect K2 at 3, 9");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 0) -
				0.9961426616053822) < eps,
			"Incorrect K2 at 4, 0");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 1) -
				0.9951625370745047) < eps,
			"Incorrect K2 at 4, 1");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 2) -
				0.9585556905472743) < eps,
			"Incorrect K2 at 4, 2");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 3) -
				0.9952971844072107) < eps,
			"Incorrect K2 at 4, 3");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 4) -
				0.9849114167202919) < eps,
			"Incorrect K2 at 4, 4");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 5) -
				0.9944417095668610) < eps,
			"Incorrect K2 at 4, 5");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 6) -
				0.9864738775087318) < eps,
			"Incorrect K2 at 4, 6");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 7) -
				0.9928713082004123) < eps,
			"Incorrect K2 at 4, 7");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 8) -
				0.9926416525703272) < eps,
			"Incorrect K2 at 4, 8");
	mu_assert(fabs(matrix_get(K2, data_2->n, data_1->n, 4, 9) -
				0.9813903448088843) < eps,
			"Incorrect K2 at 4, 9");

	free(K2);

	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data_1);
	gensvm_free_data(data_2);

	return NULL;
}

char *test_kernel_trainfactor()
{
	struct GenData *data = gensvm_init_data();

	data->n = 10;
	data->m = 5;
	data->K = 3;
	int r = 7;

	double *P = Calloc(double, data->n * r);
	double *Sigma = Calloc(double, r);

	matrix_set(P, data->n, r, 0, 0, 0.8056271362589000);
	matrix_set(P, data->n, r, 0, 1, 0.4874175854113872);
	matrix_set(P, data->n, r, 0, 2, 0.4453015882771756);
	matrix_set(P, data->n, r, 0, 3, 0.8233234072519983);
	matrix_set(P, data->n, r, 0, 4, 0.2728942849022845);
	matrix_set(P, data->n, r, 0, 5, 0.4221779040953508);
	matrix_set(P, data->n, r, 0, 6, 0.5090230605656835);
	matrix_set(P, data->n, r, 1, 0, 0.7940590105180981);
	matrix_set(P, data->n, r, 1, 1, 0.1861049005485224);
	matrix_set(P, data->n, r, 1, 2, 0.8469394287449229);
	matrix_set(P, data->n, r, 1, 3, 0.7956168453294307);
	matrix_set(P, data->n, r, 1, 4, 0.2665038412777220);
	matrix_set(P, data->n, r, 1, 5, 0.1657896146622650);
	matrix_set(P, data->n, r, 1, 6, 0.9000985899630121);
	matrix_set(P, data->n, r, 2, 0, 0.0294257611061681);
	matrix_set(P, data->n, r, 2, 1, 0.0242717976065267);
	matrix_set(P, data->n, r, 2, 2, 0.5039128672814752);
	matrix_set(P, data->n, r, 2, 3, 0.5777227081256917);
	matrix_set(P, data->n, r, 2, 4, 0.1126634767453339);
	matrix_set(P, data->n, r, 2, 5, 0.2472125660924395);
	matrix_set(P, data->n, r, 2, 6, 0.1140159498084029);
	matrix_set(P, data->n, r, 3, 0, 0.1746563833537603);
	matrix_set(P, data->n, r, 3, 1, 0.9135736087631979);
	matrix_set(P, data->n, r, 3, 2, 0.5270258081021366);
	matrix_set(P, data->n, r, 3, 3, 0.4426030703804438);
	matrix_set(P, data->n, r, 3, 4, 0.9577384085295405);
	matrix_set(P, data->n, r, 3, 5, 0.7987328321610475);
	matrix_set(P, data->n, r, 3, 6, 0.3920414766506218);
	matrix_set(P, data->n, r, 4, 0, 0.0022298761599785);
	matrix_set(P, data->n, r, 4, 1, 0.3773482059713607);
	matrix_set(P, data->n, r, 4, 2, 0.8009654729622842);
	matrix_set(P, data->n, r, 4, 3, 0.7701104553132680);
	matrix_set(P, data->n, r, 4, 4, 0.4421972387164297);
	matrix_set(P, data->n, r, 4, 5, 0.9763647809741268);
	matrix_set(P, data->n, r, 4, 6, 0.7735127328107272);
	matrix_set(P, data->n, r, 5, 0, 0.6638830667081945);
	matrix_set(P, data->n, r, 5, 1, 0.6467607601353914);
	matrix_set(P, data->n, r, 5, 2, 0.0434948735457108);
	matrix_set(P, data->n, r, 5, 3, 0.3267543833513200);
	matrix_set(P, data->n, r, 5, 4, 0.3993579061297995);
	matrix_set(P, data->n, r, 5, 5, 0.5095430570208566);
	matrix_set(P, data->n, r, 5, 6, 0.0437246682449493);
	matrix_set(P, data->n, r, 6, 0, 0.0770493004546461);
	matrix_set(P, data->n, r, 6, 1, 0.3699566427075194);
	matrix_set(P, data->n, r, 6, 2, 0.7863539761080217);
	matrix_set(P, data->n, r, 6, 3, 0.3693175185473680);
	matrix_set(P, data->n, r, 6, 4, 0.9137371110726166);
	matrix_set(P, data->n, r, 6, 5, 0.1190008864254177);
	matrix_set(P, data->n, r, 6, 6, 0.3534695686769205);
	matrix_set(P, data->n, r, 7, 0, 0.2685233952731509);
	matrix_set(P, data->n, r, 7, 1, 0.8539966432782011);
	matrix_set(P, data->n, r, 7, 2, 0.0967159557826836);
	matrix_set(P, data->n, r, 7, 3, 0.2456426390463990);
	matrix_set(P, data->n, r, 7, 4, 0.0191704319054138);
	matrix_set(P, data->n, r, 7, 5, 0.5211738041619068);
	matrix_set(P, data->n, r, 7, 6, 0.0966039750517987);
	matrix_set(P, data->n, r, 8, 0, 0.1163951898554611);
	matrix_set(P, data->n, r, 8, 1, 0.7667861436369238);
	matrix_set(P, data->n, r, 8, 2, 0.5031912600213351);
	matrix_set(P, data->n, r, 8, 3, 0.1102697774064020);
	matrix_set(P, data->n, r, 8, 4, 0.4826593547737735);
	matrix_set(P, data->n, r, 8, 5, 0.2491297759535341);
	matrix_set(P, data->n, r, 8, 6, 0.7277927065471604);
	matrix_set(P, data->n, r, 9, 0, 0.2290251898688216);
	matrix_set(P, data->n, r, 9, 1, 0.4401981048538806);
	matrix_set(P, data->n, r, 9, 2, 0.0884616753393881);
	matrix_set(P, data->n, r, 9, 3, 0.8659836346403005);
	matrix_set(P, data->n, r, 9, 4, 0.1028774221216107);
	matrix_set(P, data->n, r, 9, 5, 0.6189568847484389);
	matrix_set(P, data->n, r, 9, 6, 0.6893982492936992);

	Sigma[0] = 0.8917076715220298;
	Sigma[1] = 0.7448958325189663;
	Sigma[2] = 0.7292497207505059;
	Sigma[3] = 0.0388180555183835;
	Sigma[4] = 0.1843977014611545;
	Sigma[5] = 0.2271271152279254;
	Sigma[6] = 0.6985271697744179;


	// start test code //
	double eps = 1e-14;
	gensvm_kernel_trainfactor(data, P, Sigma, r);
	mu_assert(data->r == r, "Incorrect data->r");

	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 0, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 1) -
				0.7183838977883847) < eps,
			"Incorrect data->Z at 0, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 2) -
				0.3630753280693996) < eps,
			"Incorrect data->Z at 0, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 3) -
				0.3247360589008871) < eps,
			"Incorrect data->Z at 0, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 4) -
				0.0319598137322927) < eps,
			"Incorrect data->Z at 0, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 5) -
				0.0503210788778667) < eps,
			"Incorrect data->Z at 0, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 6) -
				0.0958880494701488) < eps,
			"Incorrect data->Z at 0, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 0, 7) -
				0.3555664378468590) < eps,
			"Incorrect data->Z at 0, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 1, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 1) -
				0.7080685113201802) < eps,
			"Incorrect data->Z at 1, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 2) -
				0.1386287648299510) < eps,
			"Incorrect data->Z at 1, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 3) -
				0.6176303419048280) < eps,
			"Incorrect data->Z at 1, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 4) -
				0.0308842988733590) < eps,
			"Incorrect data->Z at 1, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 5) -
				0.0491426957621803) < eps,
			"Incorrect data->Z at 1, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 6) -
				0.0376553169129896) < eps,
			"Incorrect data->Z at 1, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 1, 7) -
				0.6287433205648071) < eps,
			"Incorrect data->Z at 1, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 2, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 1) -
				0.0262391769187446) < eps,
			"Incorrect data->Z at 2, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 2) -
				0.0180799608848456) < eps,
			"Incorrect data->Z at 2, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 3) -
				0.3674783177476025) < eps,
			"Incorrect data->Z at 2, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 4) -
				0.0224260721582540) < eps,
			"Incorrect data->Z at 2, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 5) -
				0.0207748861504618) < eps,
			"Incorrect data->Z at 2, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 6) -
				0.0561486769846686) < eps,
			"Incorrect data->Z at 2, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 2, 7) -
				0.0796432387288058) < eps,
			"Incorrect data->Z at 2, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 3, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 1) -
				0.1557424369168406) < eps,
			"Incorrect data->Z at 3, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 2) -
				0.6805171738670187) < eps,
			"Incorrect data->Z at 3, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 3) -
				0.3843334233867928) < eps,
			"Incorrect data->Z at 3, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 4) -
				0.0171809905586351) < eps,
			"Incorrect data->Z at 3, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 5) -
				0.1766047611339114) < eps,
			"Incorrect data->Z at 3, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 6) -
				0.1814138840065695) < eps,
			"Incorrect data->Z at 3, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 3, 7) -
				0.2738516231189423) < eps,
			"Incorrect data->Z at 3, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 4, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 1) -
				0.0019883976783969) < eps,
			"Incorrect data->Z at 4, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 2) -
				0.2810851060365751) < eps,
			"Incorrect data->Z at 4, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 3) -
				0.5841038474885426) < eps,
			"Incorrect data->Z at 4, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 4) -
				0.0298941904096380) < eps,
			"Incorrect data->Z at 4, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 5) -
				0.0815401544117791) < eps,
			"Incorrect data->Z at 4, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 6) -
				0.2217589161127987) < eps,
			"Incorrect data->Z at 4, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 4, 7) -
				0.5403196600347527) < eps,
			"Incorrect data->Z at 4, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 5, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 1) -
				0.5919896235772685) < eps,
			"Incorrect data->Z at 5, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 2) -
				0.4817693948616519) < eps,
			"Incorrect data->Z at 5, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 3) -
				0.0317186243872882) < eps,
			"Incorrect data->Z at 5, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 4) -
				0.0126839697938067) < eps,
			"Incorrect data->Z at 5, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 5) -
				0.0736406799506745) < eps,
			"Incorrect data->Z at 5, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 6) -
				0.1157310446255655) < eps,
			"Incorrect data->Z at 5, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 5, 7) -
				0.0305428687584698) < eps,
			"Incorrect data->Z at 5, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 6, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 1) -
				0.0687054523008138) < eps,
			"Incorrect data->Z at 6, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 2) -
				0.2755791613655394) < eps,
			"Incorrect data->Z at 6, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 3) -
				0.5734484174878248) < eps,
			"Incorrect data->Z at 6, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 4) -
				0.0143361879388834) < eps,
			"Incorrect data->Z at 6, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 5) -
				0.1684910230215461) < eps,
			"Incorrect data->Z at 6, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 6) -
				0.0270283280433711) < eps,
			"Incorrect data->Z at 6, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 6, 7) -
				0.2469080974092735) < eps,
			"Incorrect data->Z at 6, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 7, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 1) -
				0.2394443715482110) < eps,
			"Incorrect data->Z at 7, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 2) -
				0.6361385405631184) < eps,
			"Incorrect data->Z at 7, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 3) -
				0.0705300837466403) < eps,
			"Incorrect data->Z at 7, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 4) -
				0.0095353696001854) < eps,
			"Incorrect data->Z at 7, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 5) -
				0.0035349835793759) < eps,
			"Incorrect data->Z at 7, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 6) -
				0.1183727026716577) < eps,
			"Incorrect data->Z at 7, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 7, 7) -
				0.0674805012818914) < eps,
			"Incorrect data->Z at 7, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 8, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 1) -
				0.1037904837223778) < eps,
			"Incorrect data->Z at 8, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 2) -
				0.5711758028284341) < eps,
			"Incorrect data->Z at 8, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 3) -
				0.3669520858546538) < eps,
			"Incorrect data->Z at 8, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 4) -
				0.0042804583413615) < eps,
			"Incorrect data->Z at 8, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 5) -
				0.0890012756090077) < eps,
			"Incorrect data->Z at 8, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 6) -
				0.0565841273297056) < eps,
			"Incorrect data->Z at 8, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 8, 7) -
				0.5083829794868514) < eps,
			"Incorrect data->Z at 8, 7");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 0) -
				1.0000000000000000) < eps,
			"Incorrect data->Z at 9, 0");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 1) -
				0.2042235187778177) < eps,
			"Incorrect data->Z at 9, 1");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 2) -
				0.3279017337884026) < eps,
			"Incorrect data->Z at 9, 2");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 3) -
				0.0645106520383707) < eps,
			"Incorrect data->Z at 9, 3");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 4) -
				0.0336158008074787) < eps,
			"Incorrect data->Z at 9, 4");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 5) -
				0.0189703601714739) < eps,
			"Incorrect data->Z at 9, 5");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 6) -
				0.1405818916833765) < eps,
			"Incorrect data->Z at 9, 6");
	mu_assert(fabs(matrix_get(data->Z, data->n, r+1, 9, 7) -
				0.4815634079265663) < eps,
			"Incorrect data->Z at 9, 7");

	// end test code //

	free(Sigma);
	free(P);
	gensvm_free_data(data);

	return NULL;
}

char *test_kernel_testfactor()
{
	struct GenData *test = gensvm_init_data();
	struct GenData *train = gensvm_init_data();

	int m = 3;
	train->n = 10;
	train->m = m;
	train->r = m;
	test->n = 5;
	test->m = m;

	train->Z = Calloc(double, train->n * (train->r + 1));
	train->Sigma = Calloc(double, train->r);
	double *K2 = Calloc(double, test->n * train->n);

	matrix_set(train->Z, train->n, train->r + 1, 0, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 0, 1, 0.8056271362589000);
	matrix_set(train->Z, train->n, train->r + 1, 0, 2, 0.4874175854113872);
	matrix_set(train->Z, train->n, train->r + 1, 0, 3, 0.4453015882771756);
	matrix_set(train->Z, train->n, train->r + 1, 1, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 1, 1, 0.7940590105180981);
	matrix_set(train->Z, train->n, train->r + 1, 1, 2, 0.1861049005485224);
	matrix_set(train->Z, train->n, train->r + 1, 1, 3, 0.8469394287449229);
	matrix_set(train->Z, train->n, train->r + 1, 2, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 2, 1, 0.0294257611061681);
	matrix_set(train->Z, train->n, train->r + 1, 2, 2, 0.0242717976065267);
	matrix_set(train->Z, train->n, train->r + 1, 2, 3, 0.5039128672814752);
	matrix_set(train->Z, train->n, train->r + 1, 3, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 3, 1, 0.1746563833537603);
	matrix_set(train->Z, train->n, train->r + 1, 3, 2, 0.9135736087631979);
	matrix_set(train->Z, train->n, train->r + 1, 3, 3, 0.5270258081021366);
	matrix_set(train->Z, train->n, train->r + 1, 4, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 4, 1, 0.0022298761599785);
	matrix_set(train->Z, train->n, train->r + 1, 4, 2, 0.3773482059713607);
	matrix_set(train->Z, train->n, train->r + 1, 4, 3, 0.8009654729622842);
	matrix_set(train->Z, train->n, train->r + 1, 5, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 5, 1, 0.6638830667081945);
	matrix_set(train->Z, train->n, train->r + 1, 5, 2, 0.6467607601353914);
	matrix_set(train->Z, train->n, train->r + 1, 5, 3, 0.0434948735457108);
	matrix_set(train->Z, train->n, train->r + 1, 6, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 6, 1, 0.0770493004546461);
	matrix_set(train->Z, train->n, train->r + 1, 6, 2, 0.3699566427075194);
	matrix_set(train->Z, train->n, train->r + 1, 6, 3, 0.7863539761080217);
	matrix_set(train->Z, train->n, train->r + 1, 7, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 7, 1, 0.2685233952731509);
	matrix_set(train->Z, train->n, train->r + 1, 7, 2, 0.8539966432782011);
	matrix_set(train->Z, train->n, train->r + 1, 7, 3, 0.0967159557826836);
	matrix_set(train->Z, train->n, train->r + 1, 8, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 8, 1, 0.1163951898554611);
	matrix_set(train->Z, train->n, train->r + 1, 8, 2, 0.7667861436369238);
	matrix_set(train->Z, train->n, train->r + 1, 8, 3, 0.5031912600213351);
	matrix_set(train->Z, train->n, train->r + 1, 9, 0, 1.0000000000000000);
	matrix_set(train->Z, train->n, train->r + 1, 9, 1, 0.2290251898688216);
	matrix_set(train->Z, train->n, train->r + 1, 9, 2, 0.4401981048538806);
	matrix_set(train->Z, train->n, train->r + 1, 9, 3, 0.0884616753393881);

	train->Sigma[0] = 0.3385406475521640;
	train->Sigma[1] = 0.6866094956337303;
	train->Sigma[2] = 0.2520805490786355;

	matrix_set(K2, test->n, train->n, 0, 0, 0.8233234072519983);
	matrix_set(K2, test->n, train->n, 0, 1, 0.3267543833513200);
	matrix_set(K2, test->n, train->n, 0, 2, 0.2728942849022845);
	matrix_set(K2, test->n, train->n, 0, 3, 0.3993579061297995);
	matrix_set(K2, test->n, train->n, 0, 4, 0.4221779040953508);
	matrix_set(K2, test->n, train->n, 0, 5, 0.5095430570208566);
	matrix_set(K2, test->n, train->n, 0, 6, 0.5090230605656835);
	matrix_set(K2, test->n, train->n, 0, 7, 0.0437246682449493);
	matrix_set(K2, test->n, train->n, 0, 8, 0.8917076715220298);
	matrix_set(K2, test->n, train->n, 0, 9, 0.2271271152279254);
	matrix_set(K2, test->n, train->n, 1, 0, 0.7956168453294307);
	matrix_set(K2, test->n, train->n, 1, 1, 0.3693175185473680);
	matrix_set(K2, test->n, train->n, 1, 2, 0.2665038412777220);
	matrix_set(K2, test->n, train->n, 1, 3, 0.9137371110726166);
	matrix_set(K2, test->n, train->n, 1, 4, 0.1657896146622650);
	matrix_set(K2, test->n, train->n, 1, 5, 0.1190008864254177);
	matrix_set(K2, test->n, train->n, 1, 6, 0.9000985899630121);
	matrix_set(K2, test->n, train->n, 1, 7, 0.3534695686769205);
	matrix_set(K2, test->n, train->n, 1, 8, 0.7448958325189663);
	matrix_set(K2, test->n, train->n, 1, 9, 0.6985271697744179);
	matrix_set(K2, test->n, train->n, 2, 0, 0.5777227081256917);
	matrix_set(K2, test->n, train->n, 2, 1, 0.2456426390463990);
	matrix_set(K2, test->n, train->n, 2, 2, 0.1126634767453339);
	matrix_set(K2, test->n, train->n, 2, 3, 0.0191704319054138);
	matrix_set(K2, test->n, train->n, 2, 4, 0.2472125660924395);
	matrix_set(K2, test->n, train->n, 2, 5, 0.5211738041619068);
	matrix_set(K2, test->n, train->n, 2, 6, 0.1140159498084029);
	matrix_set(K2, test->n, train->n, 2, 7, 0.0966039750517987);
	matrix_set(K2, test->n, train->n, 2, 8, 0.7292497207505059);
	matrix_set(K2, test->n, train->n, 2, 9, 0.3663483426181376);
	matrix_set(K2, test->n, train->n, 3, 0, 0.4426030703804438);
	matrix_set(K2, test->n, train->n, 3, 1, 0.1102697774064020);
	matrix_set(K2, test->n, train->n, 3, 2, 0.9577384085295405);
	matrix_set(K2, test->n, train->n, 3, 3, 0.4826593547737735);
	matrix_set(K2, test->n, train->n, 3, 4, 0.7987328321610475);
	matrix_set(K2, test->n, train->n, 3, 5, 0.2491297759535341);
	matrix_set(K2, test->n, train->n, 3, 6, 0.3920414766506218);
	matrix_set(K2, test->n, train->n, 3, 7, 0.7277927065471604);
	matrix_set(K2, test->n, train->n, 3, 8, 0.0388180555183835);
	matrix_set(K2, test->n, train->n, 3, 9, 0.0724688230975385);
	matrix_set(K2, test->n, train->n, 4, 0, 0.7701104553132680);
	matrix_set(K2, test->n, train->n, 4, 1, 0.8659836346403005);
	matrix_set(K2, test->n, train->n, 4, 2, 0.4421972387164297);
	matrix_set(K2, test->n, train->n, 4, 3, 0.1028774221216107);
	matrix_set(K2, test->n, train->n, 4, 4, 0.9763647809741268);
	matrix_set(K2, test->n, train->n, 4, 5, 0.6189568847484389);
	matrix_set(K2, test->n, train->n, 4, 6, 0.7735127328107272);
	matrix_set(K2, test->n, train->n, 4, 7, 0.6893982492936992);
	matrix_set(K2, test->n, train->n, 4, 8, 0.1843977014611545);
	matrix_set(K2, test->n, train->n, 4, 9, 0.5028223670359105);

	// start test code //
	double eps = 1e-14;
	gensvm_kernel_testfactor(test, train, K2);
	mu_assert(test->r == train->r, "Incorrect test->r");

	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 0, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 1) -
				13.4938074073234233) < eps,
			"Incorrect test->Z at 0, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 2) -
				4.9462576505461318) < eps,
			"Incorrect test->Z at 0, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 0, 3) -
				35.0141524482092095) < eps,
			"Incorrect test->Z at 0, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 1, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 1) -
				13.8904764364627322) < eps,
			"Incorrect test->Z at 1, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 2) -
				6.2592510920597659) < eps,
			"Incorrect test->Z at 1, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 1, 3) -
				40.9083467283112441) < eps,
			"Incorrect test->Z at 1, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 2, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 1) -
				10.6204422184756044) < eps,
			"Incorrect test->Z at 2, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 2) -
				3.4427837480556240) < eps,
			"Incorrect test->Z at 2, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 2, 3) -
				19.6903927993336545) < eps,
			"Incorrect test->Z at 2, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 3, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 1) -
				8.4682213303174958) < eps,
			"Incorrect test->Z at 3, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 2) -
				4.2237417955914749) < eps,
			"Incorrect test->Z at 3, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 3, 3) -
				32.7751478586422280) < eps,
			"Incorrect test->Z at 3, 3");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 0) -
				1.0000000000000000) < eps,
			"Incorrect test->Z at 4, 0");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 1) -
				18.6151364859121280) < eps,
			"Incorrect test->Z at 4, 1");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 2) -
				5.6161575228997469) < eps,
			"Incorrect test->Z at 4, 2");
	mu_assert(fabs(matrix_get(test->Z, test->n, test->r+1, 4, 3) -
				46.8107321811939130) < eps,
			"Incorrect test->Z at 4, 3");

	// end test code //

	free(K2);

	gensvm_free_data(test);
	gensvm_free_data(train);

	return NULL;
}

char *test_dsyevx()
{
	mu_test_missing();

	return NULL;
}

char *test_dlamch()
{
	mu_test_missing();

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_dot_rbf);
	mu_run_test(test_dot_poly);
	mu_run_test(test_dot_sigmoid);

	mu_run_test(test_kernel_preprocess_linear);
	mu_run_test(test_kernel_preprocess_kernel);

	mu_run_test(test_kernel_copy_kernelparam_to_data_linear);
	mu_run_test(test_kernel_copy_kernelparam_to_data_rbf);
	mu_run_test(test_kernel_copy_kernelparam_to_data_poly);
	mu_run_test(test_kernel_copy_kernelparam_to_data_sigmoid);

	mu_run_test(test_kernel_compute_rbf);
	mu_run_test(test_kernel_compute_poly);
	mu_run_test(test_kernel_compute_sigmoid);

	mu_run_test(test_kernel_eigendecomp);

	mu_run_test(test_kernel_cross_rbf);
	mu_run_test(test_kernel_cross_poly);
	mu_run_test(test_kernel_cross_sigmoid);

	mu_run_test(test_kernel_trainfactor);

	mu_run_test(test_kernel_testfactor);

	mu_run_test(test_kernel_postprocess_linear);
	mu_run_test(test_kernel_postprocess_kernel);

	return NULL;
}

RUN_TESTS(all_tests);
