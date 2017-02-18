/**
 * @file test_gensvm_train.c
 * @author G.J.J. van den Burg
 * @date 2016-12-06
 * @brief Unit tests for gensvm_train.c functions
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
#include "gensvm_train.h"

char *test_gensvm_train_seed_linear()
{
	struct GenModel *model = gensvm_init_model();
	struct GenModel *seed = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	// note that model n, m, K are set by the train function
	model->p = 1.2143;
	model->kappa = 0.90298;
	model->lambda = 0.00219038;
	model->epsilon = 1e-15;
	model->weight_idx = 1;
	model->kerneltype = K_LINEAR;

	data->n = 10;
	data->m = 3;
	data->K = 4;
	data->RAW = Calloc(double, data->n * (data->m+1));
	data->Z = data->RAW;
	data->y = Calloc(long, data->n);

	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.8056271362589000);
	matrix_set(data->Z, data->m+1, 0, 2, 0.4874175854113872);
	matrix_set(data->Z, data->m+1, 0, 3, 0.4453015882771756);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.7940590105180981);
	matrix_set(data->Z, data->m+1, 1, 2, 0.1861049005485224);
	matrix_set(data->Z, data->m+1, 1, 3, 0.8469394287449229);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, 0.0294257611061681);
	matrix_set(data->Z, data->m+1, 2, 2, 0.0242717976065267);
	matrix_set(data->Z, data->m+1, 2, 3, 0.5039128672814752);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, 0.1746563833537603);
	matrix_set(data->Z, data->m+1, 3, 2, 0.9135736087631979);
	matrix_set(data->Z, data->m+1, 3, 3, 0.5270258081021366);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, 0.0022298761599785);
	matrix_set(data->Z, data->m+1, 4, 2, 0.3773482059713607);
	matrix_set(data->Z, data->m+1, 4, 3, 0.8009654729622842);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, 0.6638830667081945);
	matrix_set(data->Z, data->m+1, 5, 2, 0.6467607601353914);
	matrix_set(data->Z, data->m+1, 5, 3, 0.0434948735457108);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, 0.0770493004546461);
	matrix_set(data->Z, data->m+1, 6, 2, 0.3699566427075194);
	matrix_set(data->Z, data->m+1, 6, 3, 0.7863539761080217);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, 0.2685233952731509);
	matrix_set(data->Z, data->m+1, 7, 2, 0.8539966432782011);
	matrix_set(data->Z, data->m+1, 7, 3, 0.0967159557826836);
	matrix_set(data->Z, data->m+1, 8, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 8, 1, 0.1163951898554611);
	matrix_set(data->Z, data->m+1, 8, 2, 0.7667861436369238);
	matrix_set(data->Z, data->m+1, 8, 3, 0.5031912600213351);
	matrix_set(data->Z, data->m+1, 9, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 9, 1, 0.2290251898688216);
	matrix_set(data->Z, data->m+1, 9, 2, 0.4401981048538806);
	matrix_set(data->Z, data->m+1, 9, 3, 0.0884616753393881);

	matrix_set(data->y, 1, 0, 0, 2);
	matrix_set(data->y, 1, 0, 1, 1);
	matrix_set(data->y, 1, 0, 2, 3);
	matrix_set(data->y, 1, 0, 3, 2);
	matrix_set(data->y, 1, 0, 4, 3);
	matrix_set(data->y, 1, 0, 5, 2);
	matrix_set(data->y, 1, 0, 6, 4);
	matrix_set(data->y, 1, 0, 7, 1);
	matrix_set(data->y, 1, 0, 8, 3);
	matrix_set(data->y, 1, 0, 9, 4);

	seed->V = Calloc(double, (data->m+1)*(data->K-1));
	matrix_set(seed->V, data->K-1, 0, 0, 0.8233234072519983);
	matrix_set(seed->V, data->K-1, 0, 1, 0.7701104553132680);
	matrix_set(seed->V, data->K-1, 0, 2, 0.1102697774064020);
	matrix_set(seed->V, data->K-1, 1, 0, 0.7956168453294307);
	matrix_set(seed->V, data->K-1, 1, 1, 0.3267543833513200);
	matrix_set(seed->V, data->K-1, 1, 2, 0.8659836346403005);
	matrix_set(seed->V, data->K-1, 2, 0, 0.5777227081256917);
	matrix_set(seed->V, data->K-1, 2, 1, 0.3693175185473680);
	matrix_set(seed->V, data->K-1, 2, 2, 0.2728942849022845);
	matrix_set(seed->V, data->K-1, 3, 0, 0.4426030703804438);
	matrix_set(seed->V, data->K-1, 3, 1, 0.2456426390463990);
	matrix_set(seed->V, data->K-1, 3, 2, 0.2665038412777220);

	// start test code //
	gensvm_train(model, data, seed);

	mu_assert(model->n == data->n, "Incorrect model n");
	mu_assert(model->m == data->m, "Incorrect model m");
	mu_assert(model->K == data->K, "Incorrect model K");

	double eps = 1e-13;
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 0) -
				-1.1907736868272805) < eps,
			"Incorrect model->V at 0, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 1) -
				1.8651287814979396) < eps,
			"Incorrect model->V at 0, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 0, 2) -
				1.7250030581662932) < eps,
			"Incorrect model->V at 0, 2");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 0) -
				0.7925100058806183) < eps,
			"Incorrect model->V at 1, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 1) -
				-3.6093428916761665) < eps,
			"Incorrect model->V at 1, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 1, 2) -
				-1.3394018960329377) < eps,
			"Incorrect model->V at 1, 2");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 0) -
				1.5203132433193016) < eps,
			"Incorrect model->V at 2, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 1) -
				-1.9118604362643852) < eps,
			"Incorrect model->V at 2, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 2, 2) -
				-1.7939246097629342) < eps,
			"Incorrect model->V at 2, 2");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 0) -
				0.0658817457370326) < eps,
			"Incorrect model->V at 3, 0");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 1) -
				0.6547924025329720) < eps,
			"Incorrect model->V at 3, 1");
	mu_assert(fabs(matrix_get(model->V, model->K-1, 3, 2) -
				-0.6773346708737853) < eps,
			"Incorrect model->V at 3, 2");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_model(seed);
	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_train_seed_kernel()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	// note that model n, m, K are set by the train function
	model->p = 1.2143;
	model->kappa = 0.90298;
	model->lambda = 0.00219038;
	model->epsilon = 1e-15;
	model->weight_idx = 1;
	model->kerneltype = K_RBF;
	model->gamma = 0.348;
	model->kernel_eigen_cutoff = 5e-3;

	data->n = 10;
	data->m = 5;
	data->K = 4;
	data->RAW = Calloc(double, data->n * (data->m+1));
	data->Z = data->RAW;
	matrix_set(data->Z, data->m+1, 0, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 0, 1, 0.0657799204744603);
	matrix_set(data->Z, data->m+1, 0, 2, 0.2576653302581353);
	matrix_set(data->Z, data->m+1, 0, 3, 0.0221000752651170);
	matrix_set(data->Z, data->m+1, 0, 4, 0.6666929354133441);
	matrix_set(data->Z, data->m+1, 0, 5, 0.6178892590244618);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 1, 1, 0.9797668012781366);
	matrix_set(data->Z, data->m+1, 1, 2, 0.7636361573939686);
	matrix_set(data->Z, data->m+1, 1, 3, 0.3195806959299131);
	matrix_set(data->Z, data->m+1, 1, 4, 0.2947771273705799);
	matrix_set(data->Z, data->m+1, 1, 5, 0.8358899802514324);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 2, 1, 0.9473849700145257);
	matrix_set(data->Z, data->m+1, 2, 2, 0.8682867844262768);
	matrix_set(data->Z, data->m+1, 2, 3, 0.7116177283612393);
	matrix_set(data->Z, data->m+1, 2, 4, 0.5092752476335579);
	matrix_set(data->Z, data->m+1, 2, 5, 0.1046097156193449);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 3, 1, 0.5846585351601830);
	matrix_set(data->Z, data->m+1, 3, 2, 0.4076887966131124);
	matrix_set(data->Z, data->m+1, 3, 3, 0.8661556045821296);
	matrix_set(data->Z, data->m+1, 3, 4, 0.0904082115920005);
	matrix_set(data->Z, data->m+1, 3, 5, 0.0799888711622944);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 4, 1, 0.8112201081242789);
	matrix_set(data->Z, data->m+1, 4, 2, 0.3112642417912803);
	matrix_set(data->Z, data->m+1, 4, 3, 0.7902557587124555);
	matrix_set(data->Z, data->m+1, 4, 4, 0.3001992968661185);
	matrix_set(data->Z, data->m+1, 4, 5, 0.6030590437920392);
	matrix_set(data->Z, data->m+1, 5, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 5, 1, 0.0098576324913424);
	matrix_set(data->Z, data->m+1, 5, 2, 0.5686603332895077);
	matrix_set(data->Z, data->m+1, 5, 3, 0.9933970661175713);
	matrix_set(data->Z, data->m+1, 5, 4, 0.5215400841900655);
	matrix_set(data->Z, data->m+1, 5, 5, 0.4307310515440625);
	matrix_set(data->Z, data->m+1, 6, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 6, 1, 0.2773296707204919);
	matrix_set(data->Z, data->m+1, 6, 2, 0.5114254316901164);
	matrix_set(data->Z, data->m+1, 6, 3, 0.5057613745592034);
	matrix_set(data->Z, data->m+1, 6, 4, 0.6411421568717217);
	matrix_set(data->Z, data->m+1, 6, 5, 0.3114658800558432);
	matrix_set(data->Z, data->m+1, 7, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 7, 1, 0.7195909422652624);
	matrix_set(data->Z, data->m+1, 7, 2, 0.7754155342547566);
	matrix_set(data->Z, data->m+1, 7, 3, 0.5955643008534165);
	matrix_set(data->Z, data->m+1, 7, 4, 0.5920949759391909);
	matrix_set(data->Z, data->m+1, 7, 5, 0.7029537245575100);
	matrix_set(data->Z, data->m+1, 8, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 8, 1, 0.3792168380438625);
	matrix_set(data->Z, data->m+1, 8, 2, 0.1920178667928286);
	matrix_set(data->Z, data->m+1, 8, 3, 0.2742847467912714);
	matrix_set(data->Z, data->m+1, 8, 4, 0.2337979820454409);
	matrix_set(data->Z, data->m+1, 8, 5, 0.3978991644742557);
	matrix_set(data->Z, data->m+1, 9, 0, 1.0000000000000000);
	matrix_set(data->Z, data->m+1, 9, 1, 0.0797813938980598);
	matrix_set(data->Z, data->m+1, 9, 2, 0.5863311792537960);
	matrix_set(data->Z, data->m+1, 9, 3, 0.8565105304166337);
	matrix_set(data->Z, data->m+1, 9, 4, 0.8266471128109379);
	matrix_set(data->Z, data->m+1, 9, 5, 0.8070610088865674);

	data->y = Calloc(long, data->n);
	matrix_set(data->y, 1, 0, 0, 2);
	matrix_set(data->y, 1, 0, 1, 1);
	matrix_set(data->y, 1, 0, 2, 3);
	matrix_set(data->y, 1, 0, 3, 2);
	matrix_set(data->y, 1, 0, 4, 3);
	matrix_set(data->y, 1, 0, 5, 2);
	matrix_set(data->y, 1, 0, 6, 4);
	matrix_set(data->y, 1, 0, 7, 1);
	matrix_set(data->y, 1, 0, 8, 3);
	matrix_set(data->y, 1, 0, 9, 4);

	// start test code //

	// because the kernel eigendecomposition isn't known in advance, 
	// there's no way to seed the model when using kernels. We therefore 
	// use seed == NULL here. Note that due to the Memset in 
	// gensvm_reallocate_model(), V will be a matrix of zeros after 
	// reallocation, so we compare with the V = 0 result from Octave.
	gensvm_train(model, data, NULL);

	mu_assert(model->n == data->n, "Incorrect model n");
	mu_assert(model->m == data->r, "Incorrect model m");
	mu_assert(model->K == data->K, "Incorrect model K");

	double eps = 1e-13;
	mu_assert(fabs(matrix_get(data->Sigma, 1, 0, 0) -
				7.8302939172918506) < eps,
			"Incorrect data->Sigma at 0, 0");
	mu_assert(fabs(matrix_get(data->Sigma, 1, 1, 0) -
				0.7947913383766066) < eps,
			"Incorrect data->Sigma at 1, 0");
	mu_assert(fabs(matrix_get(data->Sigma, 1, 2, 0) -
				0.5288740088908547) < eps,
			"Incorrect data->Sigma at 2, 0");
	mu_assert(fabs(matrix_get(data->Sigma, 1, 3, 0) -
				0.4537982052555444) < eps,
			"Incorrect data->Sigma at 3, 0");
	mu_assert(fabs(matrix_get(data->Sigma, 1, 4, 0) -
				0.2226012271232192) < eps,
			"Incorrect data->Sigma at 4, 0");
	mu_assert(fabs(matrix_get(data->Sigma, 1, 5, 0) -
				0.0743004417495061) < eps,
			"Incorrect data->Sigma at 5, 0");

	// we need a large eps here because there are numerical precision 
	// differences between the C and Octave implementations. We also 
	// compare with absolute values because of variability in the 
	// eigendecomposition.
	eps = 1e-7;
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 0, 0)) -
				fabs(5.0555413160638665)) < eps,
			"Incorrect model->V at 0, 0");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 0, 1)) -
				fabs(-2.2586632211763198)) < eps,
			"Incorrect model->V at 0, 1");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 0, 2)) -
				fabs(-4.5572671806963143)) < eps,
			"Incorrect model->V at 0, 2");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 1, 0)) -
				fabs(-1.9627432869558412)) < eps,
			"Incorrect model->V at 1, 0");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 1, 1)) -
				fabs(0.9934555242449399)) < eps,
			"Incorrect model->V at 1, 1");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 1, 2)) -
				fabs(1.7855287218670219)) < eps,
			"Incorrect model->V at 1, 2");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 2, 0)) -
				fabs(1.9393083227054353)) < eps,
			"Incorrect model->V at 2, 0");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 2, 1)) -
				fabs(-1.1958487809502740)) < eps,
			"Incorrect model->V at 2, 1");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 2, 2)) -
				fabs(2.1140967864804359)) < eps,
			"Incorrect model->V at 2, 2");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 3, 0)) -
				fabs(2.3909204618652535)) < eps,
			"Incorrect model->V at 3, 0");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 3, 1)) -
				fabs(-0.2834554569573399)) < eps,
			"Incorrect model->V at 3, 1");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 3, 2)) -
				fabs(1.0926232371314393)) < eps,
			"Incorrect model->V at 3, 2");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 4, 0)) -
				fabs(3.3374545494113272)) < eps,
			"Incorrect model->V at 4, 0");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 4, 1)) -
				fabs(1.6699291195221897)) < eps,
			"Incorrect model->V at 4, 1");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 4, 2)) -
				fabs(-1.4345249893609275)) < eps,
			"Incorrect model->V at 4, 2");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 5, 0)) -
				fabs(-0.0221825925355533)) < eps,
			"Incorrect model->V at 5, 0");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 5, 1)) -
				fabs(-0.1216077739550210)) < eps,
			"Incorrect model->V at 5, 1");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 5, 2)) -
				fabs(-0.7900947982642630)) < eps,
			"Incorrect model->V at 5, 2");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 6, 0)) -
				fabs(-0.0076471781062262)) < eps,
			"Incorrect model->V at 6, 0");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 6, 1)) -
				fabs(-0.8781872510019056)) < eps,
			"Incorrect model->V at 6, 1");
	mu_assert(fabs(fabs(matrix_get(model->V, model->K-1, 6, 2)) -
				fabs(-0.2782284589344380)) < eps,
			"Incorrect model->V at 6, 2");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();

	mu_run_test(test_gensvm_train_seed_linear);
	mu_run_test(test_gensvm_train_seed_kernel);

	return NULL;
}

RUN_TESTS(all_tests);
