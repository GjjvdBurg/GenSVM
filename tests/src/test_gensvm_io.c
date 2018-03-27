/**
 * @file test_gensvm_io.c
 * @author G.J.J. van den Burg
 * @date 2016-09-01
 * @brief Unit tests for gensvm_io.c functions
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
#include "gensvm_io.h"

char *test_gensvm_read_data()
{
	char *filename = "./data/test_file_read_data.txt";
	struct GenData *data = gensvm_init_data();

	// start test code //
	gensvm_read_data(data, filename);

	// check if dimensions are correctly read
	mu_assert(data->n == 5, "Incorrect value for n");
	mu_assert(data->m == 3, "Incorrect value for m");
	mu_assert(data->r == 3, "Incorrect value for r");
	mu_assert(data->K == 4, "Incorrect value for K");

	// check if all data is read correctly.
	mu_assert(matrix_get(data->Z, data->m+1, 0, 1) == 0.7065937536993949,
			"Incorrect Z value at 0, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 2) == 0.7016517970438980,
			"Incorrect Z value at 0, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 3) == 0.1548611397288129,
			"Incorrect Z value at 0, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 1) == 0.4604987687863951,
			"Incorrect Z value at 1, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 2) == 0.6374142980176117,
			"Incorrect Z value at 1, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 3) == 0.0370930278245423,
			"Incorrect Z value at 1, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 1) == 0.3798777132278375,
			"Incorrect Z value at 2, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 2) == 0.5745070018747664,
			"Incorrect Z value at 2, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 3) == 0.2570906697837264,
			"Incorrect Z value at 2, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 1) == 0.2789376050039792,
			"Incorrect Z value at 3, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 2) == 0.4853242744610165,
			"Incorrect Z value at 3, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 3) == 0.1894010436762711,
			"Incorrect Z value at 3, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 1) == 0.7630904372339489,
			"Incorrect Z value at 4, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 2) == 0.1341546320318005,
			"Incorrect Z value at 4, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 3) == 0.6827430912944857,
			"Incorrect Z value at 4, 3");
	// check if RAW = Z
	mu_assert(data->Z == data->RAW, "Z pointer doesn't equal RAW pointer");

	// check if labels read correctly
	mu_assert(data->y[0] == 2, "Incorrect label read at 0");
	mu_assert(data->y[1] == 1, "Incorrect label read at 1");
	mu_assert(data->y[2] == 3, "Incorrect label read at 2");
	mu_assert(data->y[3] == 4, "Incorrect label read at 3");
	mu_assert(data->y[4] == 3, "Incorrect label read at 4");

	// check if the column of ones is added
	mu_assert(matrix_get(data->Z, data->m+1, 0, 0) == 1,
			"Incorrect Z value at 0, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 0) == 1,
			"Incorrect Z value at 1, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 0) == 1,
			"Incorrect Z value at 2, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 0) == 1,
			"Incorrect Z value at 3, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 0) == 1,
			"Incorrect Z value at 4, 0");

	// end test code //

	gensvm_free_data(data);
	return NULL;
}

char *test_gensvm_read_data_sparse()
{
	char *filename = "./data/test_file_read_data_sparse.txt";
	struct GenData *data = gensvm_init_data();

	// start test code //
	gensvm_read_data(data, filename);

	// check if dimensions are correctly read
	mu_assert(data->n == 10, "Incorrect value for n");
	mu_assert(data->m == 3, "Incorrect value for m");
	mu_assert(data->r == 3, "Incorrect value for r");
	mu_assert(data->K == 4, "Incorrect value for K");

	// check if dense data pointers are NULL
	mu_assert(data->Z == NULL, "Z pointer isn't NULL");
	mu_assert(data->RAW == NULL, "RAW pointer isn't NULL");

	// check sparse data structure
	mu_assert(data->spZ != NULL, "spZ is NULL");
	mu_assert(data->spZ->nnz == 14, "Incorrect nnz");
	mu_assert(data->spZ->n_row == 10, "Incorrect n_row");
	mu_assert(data->spZ->n_col == 4, "Incorrect n_col");

	// check sparse values
	mu_assert(data->spZ->values[0] == 1.0,
			"Incorrect nonzero value at 0");
	mu_assert(data->spZ->values[1] == 0.7016517970438980,
			"Incorrect nonzero value at 1");
	mu_assert(data->spZ->values[2] == 1.0,
			"Incorrect nonzero value at 2");
	mu_assert(data->spZ->values[3] == 0.0370930278245423,
			"Incorrect nonzero value at 3");
	mu_assert(data->spZ->values[4] == 1.0,
			"Incorrect nonzero value at 4");
	mu_assert(data->spZ->values[5] == 1.0,
			"Incorrect nonzero value at 5");
	mu_assert(data->spZ->values[6] == 0.4853242744610165,
			"Incorrect nonzero value at 6");
	mu_assert(data->spZ->values[7] == 1.0,
			"Incorrect nonzero value at 7");
	mu_assert(data->spZ->values[8] == 0.7630904372339489,
			"Incorrect nonzero value at 8");
	mu_assert(data->spZ->values[9] == 1.0,
			"Incorrect nonzero value at 9");
	mu_assert(data->spZ->values[10] == 1.0,
			"Incorrect nonzero value at 10");
	mu_assert(data->spZ->values[11] == 1.0,
			"Incorrect nonzero value at 11");
	mu_assert(data->spZ->values[12] == 1.0,
			"Incorrect nonzero value at 12");
	mu_assert(data->spZ->values[13] == 1.0,
			"Incorrect nonzero value at 13");

	// check sparse row lengths
	mu_assert(data->spZ->ia[0] == 0, "Incorrect ia value at 0");
	mu_assert(data->spZ->ia[1] == 2, "Incorrect ia value at 1");
	mu_assert(data->spZ->ia[2] == 4, "Incorrect ia value at 2");
	mu_assert(data->spZ->ia[3] == 5, "Incorrect ia value at 3");
	mu_assert(data->spZ->ia[4] == 7, "Incorrect ia value at 4");
	mu_assert(data->spZ->ia[5] == 9, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[6] == 10, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[7] == 11, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[8] == 12, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[9] == 13, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[10] == 14, "Incorrect ia value at 5");

	// check sparse column indices
	mu_assert(data->spZ->ja[0] == 0, "Incorrect ja value at 0");
	mu_assert(data->spZ->ja[1] == 2, "Incorrect ja value at 1");
	mu_assert(data->spZ->ja[2] == 0, "Incorrect ja value at 2");
	mu_assert(data->spZ->ja[3] == 3, "Incorrect ja value at 3");
	mu_assert(data->spZ->ja[4] == 0, "Incorrect ja value at 4");
	mu_assert(data->spZ->ja[5] == 0, "Incorrect ja value at 5");
	mu_assert(data->spZ->ja[6] == 2, "Incorrect ja value at 6");
	mu_assert(data->spZ->ja[7] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[8] == 1, "Incorrect ja value at 8");
	mu_assert(data->spZ->ja[9] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[10] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[11] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[12] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[13] == 0, "Incorrect ja value at 7");

	// check if labels read correctly
	mu_assert(data->y[0] == 2, "Incorrect label read at 0");
	mu_assert(data->y[1] == 1, "Incorrect label read at 1");
	mu_assert(data->y[2] == 3, "Incorrect label read at 2");
	mu_assert(data->y[3] == 4, "Incorrect label read at 3");
	mu_assert(data->y[4] == 3, "Incorrect label read at 4");

	// end test code //

	gensvm_free_data(data);
	return NULL;
}

char *test_gensvm_read_data_no_label()
{
	char *filename = "./data/test_file_read_data_no_label.txt";
	struct GenData *data = gensvm_init_data();

	// start test code //
	gensvm_read_data(data, filename);

	// check if dimensions are correctly read
	mu_assert(data->n == 5, "Incorrect value for n");
	mu_assert(data->m == 3, "Incorrect value for m");
	mu_assert(data->r == 3, "Incorrect value for r");
	mu_assert(data->K == 0, "Incorrect value for K");

	// check if all data is read correctly.
	mu_assert(matrix_get(data->Z, data->m+1, 0, 1) == 0.7065937536993949,
			"Incorrect Z value at 0, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 2) == 0.7016517970438980,
			"Incorrect Z value at 0, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 3) == 0.1548611397288129,
			"Incorrect Z value at 0, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 1) == 0.4604987687863951,
			"Incorrect Z value at 1, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 2) == 0.6374142980176117,
			"Incorrect Z value at 1, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 3) == 0.0370930278245423,
			"Incorrect Z value at 1, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 1) == 0.3798777132278375,
			"Incorrect Z value at 2, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 2) == 0.5745070018747664,
			"Incorrect Z value at 2, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 3) == 0.2570906697837264,
			"Incorrect Z value at 2, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 1) == 0.2789376050039792,
			"Incorrect Z value at 3, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 2) == 0.4853242744610165,
			"Incorrect Z value at 3, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 3) == 0.1894010436762711,
			"Incorrect Z value at 3, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 1) == 0.7630904372339489,
			"Incorrect Z value at 4, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 2) == 0.1341546320318005,
			"Incorrect Z value at 4, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 3) == 0.6827430912944857,
			"Incorrect Z value at 4, 3");
	// check if RAW = Z
	mu_assert(data->Z == data->RAW, "Z pointer doesn't equal RAW pointer");

	// check if labels read correctly
	mu_assert(data->y == NULL, "Outcome pointer is not NULL");

	// check if the column of ones is added
	mu_assert(matrix_get(data->Z, data->m+1, 0, 0) == 1,
			"Incorrect Z value at 0, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 0) == 1,
			"Incorrect Z value at 1, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 0) == 1,
			"Incorrect Z value at 2, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 0) == 1,
			"Incorrect Z value at 3, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 0) == 1,
			"Incorrect Z value at 4, 0");

	// end test code //

	gensvm_free_data(data);
	return NULL;
}

char *test_gensvm_read_data_libsvm()
{
	char *filename = "./data/test_file_read_data_libsvm.txt";
	struct GenData *data = gensvm_init_data();

	// start test code //
	gensvm_read_data_libsvm(data, filename);

	// check if dimensions are correctly read
	mu_assert(data->n == 5, "Incorrect value for n");
	mu_assert(data->m == 3, "Incorrect value for m");
	mu_assert(data->r == 3, "Incorrect value for r");
	mu_assert(data->K == 4, "Incorrect value for K");

	// check if all data is read correctly.
	mu_assert(matrix_get(data->Z, data->m+1, 0, 1) == 0.7065937536993949,
			"Incorrect Z value at 0, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 2) == 0.7016517970438980,
			"Incorrect Z value at 0, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 3) == 0.1548611397288129,
			"Incorrect Z value at 0, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 1) == 0.4604987687863951,
			"Incorrect Z value at 1, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 2) == 0.6374142980176117,
			"Incorrect Z value at 1, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 3) == 0.0370930278245423,
			"Incorrect Z value at 1, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 1) == 0.3798777132278375,
			"Incorrect Z value at 2, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 2) == 0.5745070018747664,
			"Incorrect Z value at 2, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 3) == 0.2570906697837264,
			"Incorrect Z value at 2, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 1) == 0.2789376050039792,
			"Incorrect Z value at 3, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 2) == 0.4853242744610165,
			"Incorrect Z value at 3, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 3) == 0.1894010436762711,
			"Incorrect Z value at 3, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 1) == 0.7630904372339489,
			"Incorrect Z value at 4, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 2) == 0.1341546320318005,
			"Incorrect Z value at 4, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 3) == 0.6827430912944857,
			"Incorrect Z value at 4, 3");
	// check if RAW = Z
	mu_assert(data->Z == data->RAW, "Z pointer doesn't equal RAW pointer");

	// check if labels read correctly
	mu_assert(data->y[0] == 2, "Incorrect label read at 0");
	mu_assert(data->y[1] == 1, "Incorrect label read at 1");
	mu_assert(data->y[2] == 3, "Incorrect label read at 2");
	mu_assert(data->y[3] == 4, "Incorrect label read at 3");
	mu_assert(data->y[4] == 3, "Incorrect label read at 4");

	// check if the column of ones is added
	mu_assert(matrix_get(data->Z, data->m+1, 0, 0) == 1,
			"Incorrect Z value at 0, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 0) == 1,
			"Incorrect Z value at 1, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 0) == 1,
			"Incorrect Z value at 2, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 0) == 1,
			"Incorrect Z value at 3, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 0) == 1,
			"Incorrect Z value at 4, 0");

	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_read_data_libsvm_0based()
{
	char *filename = "./data/test_file_read_data_libsvm_0.txt";
	struct GenData *data = gensvm_init_data();

	// start test code //
	gensvm_read_data_libsvm(data, filename);

	// check if dimensions are correctly read
	mu_assert(data->n == 5, "Incorrect value for n");
	mu_assert(data->m == 3, "Incorrect value for m");
	mu_assert(data->r == 3, "Incorrect value for r");
	mu_assert(data->K == 4, "Incorrect value for K");

	// check if all data is read correctly.
	mu_assert(matrix_get(data->Z, data->m+1, 0, 1) == 0.7065937536993949,
			"Incorrect Z value at 0, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 2) == 0.7016517970438980,
			"Incorrect Z value at 0, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 3) == 0.1548611397288129,
			"Incorrect Z value at 0, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 1) == 0.4604987687863951,
			"Incorrect Z value at 1, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 2) == 0.6374142980176117,
			"Incorrect Z value at 1, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 3) == 0.0370930278245423,
			"Incorrect Z value at 1, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 1) == 0.3798777132278375,
			"Incorrect Z value at 2, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 2) == 0.5745070018747664,
			"Incorrect Z value at 2, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 3) == 0.2570906697837264,
			"Incorrect Z value at 2, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 1) == 0.2789376050039792,
			"Incorrect Z value at 3, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 2) == 0.4853242744610165,
			"Incorrect Z value at 3, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 3) == 0.1894010436762711,
			"Incorrect Z value at 3, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 1) == 0.7630904372339489,
			"Incorrect Z value at 4, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 2) == 0.1341546320318005,
			"Incorrect Z value at 4, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 3) == 0.6827430912944857,
			"Incorrect Z value at 4, 3");
	// check if RAW = Z
	mu_assert(data->Z == data->RAW, "Z pointer doesn't equal RAW pointer");

	// check if labels read correctly
	mu_assert(data->y[0] == 2, "Incorrect label read at 0");
	mu_assert(data->y[1] == 1, "Incorrect label read at 1");
	mu_assert(data->y[2] == 3, "Incorrect label read at 2");
	mu_assert(data->y[3] == 4, "Incorrect label read at 3");
	mu_assert(data->y[4] == 3, "Incorrect label read at 4");

	// check if the column of ones is added
	mu_assert(matrix_get(data->Z, data->m+1, 0, 0) == 1,
			"Incorrect Z value at 0, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 0) == 1,
			"Incorrect Z value at 1, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 0) == 1,
			"Incorrect Z value at 2, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 0) == 1,
			"Incorrect Z value at 3, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 0) == 1,
			"Incorrect Z value at 4, 0");

	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_read_data_libsvm_sparse()
{
	char *filename = "./data/test_file_read_data_sparse_libsvm.txt";
	struct GenData *data = gensvm_init_data();

	// start test code //
	gensvm_read_data_libsvm(data, filename);

	// check if dimensions are correctly read
	mu_assert(data->n == 10, "Incorrect value for n");
	mu_assert(data->m == 3, "Incorrect value for m");
	mu_assert(data->r == 3, "Incorrect value for r");
	mu_assert(data->K == 4, "Incorrect value for K");

	// check if dense data pointers are NULL
	mu_assert(data->Z == NULL, "Z pointer isn't NULL");
	mu_assert(data->RAW == NULL, "RAW pointer isn't NULL");

	// check sparse data structure
	mu_assert(data->spZ != NULL, "spZ is NULL");
	mu_assert(data->spZ->nnz == 14, "Incorrect nnz");
	mu_assert(data->spZ->n_row == 10, "Incorrect n_row");
	mu_assert(data->spZ->n_col == 4, "Incorrect n_col");

	// check sparse values
	mu_assert(data->spZ->values[0] == 1.0,
			"Incorrect nonzero value at 0");
	mu_assert(data->spZ->values[1] == 0.7016517970438980,
			"Incorrect nonzero value at 1");
	mu_assert(data->spZ->values[2] == 1.0,
			"Incorrect nonzero value at 2");
	mu_assert(data->spZ->values[3] == 0.0370930278245423,
			"Incorrect nonzero value at 3");
	mu_assert(data->spZ->values[4] == 1.0,
			"Incorrect nonzero value at 4");
	mu_assert(data->spZ->values[5] == 1.0,
			"Incorrect nonzero value at 5");
	mu_assert(data->spZ->values[6] == 0.4853242744610165,
			"Incorrect nonzero value at 6");
	mu_assert(data->spZ->values[7] == 1.0,
			"Incorrect nonzero value at 7");
	mu_assert(data->spZ->values[8] == 0.7630904372339489,
			"Incorrect nonzero value at 8");
	mu_assert(data->spZ->values[9] == 1.0,
			"Incorrect nonzero value at 9");
	mu_assert(data->spZ->values[10] == 1.0,
			"Incorrect nonzero value at 10");
	mu_assert(data->spZ->values[11] == 1.0,
			"Incorrect nonzero value at 11");
	mu_assert(data->spZ->values[12] == 1.0,
			"Incorrect nonzero value at 12");
	mu_assert(data->spZ->values[13] == 1.0,
			"Incorrect nonzero value at 13");

	// check sparse row lengths
	mu_assert(data->spZ->ia[0] == 0, "Incorrect ia value at 0");
	mu_assert(data->spZ->ia[1] == 2, "Incorrect ia value at 1");
	mu_assert(data->spZ->ia[2] == 4, "Incorrect ia value at 2");
	mu_assert(data->spZ->ia[3] == 5, "Incorrect ia value at 3");
	mu_assert(data->spZ->ia[4] == 7, "Incorrect ia value at 4");
	mu_assert(data->spZ->ia[5] == 9, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[6] == 10, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[7] == 11, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[8] == 12, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[9] == 13, "Incorrect ia value at 5");
	mu_assert(data->spZ->ia[10] == 14, "Incorrect ia value at 5");

	// check sparse column indices
	mu_assert(data->spZ->ja[0] == 0, "Incorrect ja value at 0");
	mu_assert(data->spZ->ja[1] == 2, "Incorrect ja value at 1");
	mu_assert(data->spZ->ja[2] == 0, "Incorrect ja value at 2");
	mu_assert(data->spZ->ja[3] == 3, "Incorrect ja value at 3");
	mu_assert(data->spZ->ja[4] == 0, "Incorrect ja value at 4");
	mu_assert(data->spZ->ja[5] == 0, "Incorrect ja value at 5");
	mu_assert(data->spZ->ja[6] == 2, "Incorrect ja value at 6");
	mu_assert(data->spZ->ja[7] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[8] == 1, "Incorrect ja value at 8");
	mu_assert(data->spZ->ja[9] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[10] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[11] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[12] == 0, "Incorrect ja value at 7");
	mu_assert(data->spZ->ja[13] == 0, "Incorrect ja value at 7");

	// check if labels read correctly
	mu_assert(data->y[0] == 2, "Incorrect label read at 0");
	mu_assert(data->y[1] == 1, "Incorrect label read at 1");
	mu_assert(data->y[2] == 3, "Incorrect label read at 2");
	mu_assert(data->y[3] == 4, "Incorrect label read at 3");
	mu_assert(data->y[4] == 3, "Incorrect label read at 4");

	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_read_data_libsvm_no_label()
{
	char *filename = "./data/test_file_read_data_no_label_libsvm.txt";
	struct GenData *data = gensvm_init_data();

	// start test code //
	gensvm_read_data_libsvm(data, filename);

	// check if dimensions are correctly read
	mu_assert(data->n == 5, "Incorrect value for n");
	mu_assert(data->m == 3, "Incorrect value for m");
	mu_assert(data->r == 3, "Incorrect value for r");
	mu_assert(data->K == 0, "Incorrect value for K");

	// check if all data is read correctly.
	mu_assert(matrix_get(data->Z, data->m+1, 0, 1) == 0.7065937536993949,
			"Incorrect Z value at 0, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 2) == 0.7016517970438980,
			"Incorrect Z value at 0, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 0, 3) == 0.1548611397288129,
			"Incorrect Z value at 0, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 1) == 0.4604987687863951,
			"Incorrect Z value at 1, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 2) == 0.6374142980176117,
			"Incorrect Z value at 1, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 3) == 0.0370930278245423,
			"Incorrect Z value at 1, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 1) == 0.3798777132278375,
			"Incorrect Z value at 2, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 2) == 0.5745070018747664,
			"Incorrect Z value at 2, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 3) == 0.2570906697837264,
			"Incorrect Z value at 2, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 1) == 0.2789376050039792,
			"Incorrect Z value at 3, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 2) == 0.4853242744610165,
			"Incorrect Z value at 3, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 3) == 0.1894010436762711,
			"Incorrect Z value at 3, 3");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 1) == 0.7630904372339489,
			"Incorrect Z value at 4, 1");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 2) == 0.1341546320318005,
			"Incorrect Z value at 4, 2");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 3) == 0.6827430912944857,
			"Incorrect Z value at 4, 3");
	// check if RAW = Z
	mu_assert(data->Z == data->RAW, "Z pointer doesn't equal RAW pointer");

	// check if labels read correctly
	mu_assert(data->y == NULL, "Outcome pointer is not NULL");

	// check if the column of ones is added
	mu_assert(matrix_get(data->Z, data->m+1, 0, 0) == 1,
			"Incorrect Z value at 0, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 1, 0) == 1,
			"Incorrect Z value at 1, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 2, 0) == 1,
			"Incorrect Z value at 2, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 3, 0) == 1,
			"Incorrect Z value at 3, 0");
	mu_assert(matrix_get(data->Z, data->m+1, 4, 0) == 1,
			"Incorrect Z value at 4, 0");

	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *test_gensvm_read_model()
{
	struct GenModel *model = gensvm_init_model();
	char *filename = "./data/test_read_model.txt";

	// start test code //
	gensvm_read_model(model, filename);

	mu_assert(model->p == 0.21398, "Incorrect read for model->p");
	mu_assert(model->lambda == 0.902131, "Incorrect read for "
			"model->lambda");
	mu_assert(model->kappa == 1.0213, "Incorrect read for model->kappa");
	mu_assert(model->epsilon == 1e-10, "Incorrect read for "
			"model->epsilon");
	mu_assert(model->weight_idx == 2, "Incorrect read for "
			"model->weight_idx");

	mu_assert(strcmp(model->data_file, "./data/test_file_read_data.txt")
			== 0, "Incorrect read for model->data_file");
	mu_assert(model->n == 10, "Incorrect read for model->n");
	mu_assert(model->m == 2, "Incorrect read for model->m");
	mu_assert(model->K == 3, "Incorrect read for model->K");

	mu_assert(matrix_get(model->V, model->K-1, 0, 0) == 0.1234,
			"Incorrect model->V element at 0, 0");
	mu_assert(matrix_get(model->V, model->K-1, 0, 1) == -0.4321,
			"Incorrect model->V element at 0, 1");
	mu_assert(matrix_get(model->V, model->K-1, 1, 0) == -0.5678,
			"Incorrect model->V element at 1, 0");
	mu_assert(matrix_get(model->V, model->K-1, 1, 1) == 0.9876,
			"Incorrect model->V element at 1, 1");
	mu_assert(matrix_get(model->V, model->K-1, 2, 0) == 0.9012,
			"Incorrect model->V element at 2, 0");
	mu_assert(matrix_get(model->V, model->K-1, 2, 1) == -0.5555,
			"Incorrect model->V element at 2, 1");

	// end test code //

	gensvm_free_model(model);

	return NULL;
}

char *test_gensvm_write_model()
{
	struct GenModel *model = gensvm_init_model();

	model->p = 0.90328;
	model->lambda = 0.0130;
	model->kappa = 1.1832;
	model->epsilon = 1e-8;
	model->weight_idx = 1;
	model->data_file = strdup("./data/test_file_read_data.txt");
	model->n = 10;
	model->m = 2;
	model->K = 3;

	model->V = Calloc(double, (model->m+1)*(model->K-1));
	matrix_set(model->V, model->K-1, 0, 0, 0.4989893785603748);
	matrix_set(model->V, model->K-1, 0, 1, 0.0599082796573645);
	matrix_set(model->V, model->K-1, 1, 0, 0.7918204761759593);
	matrix_set(model->V, model->K-1, 1, 1, 0.6456613497110559);
	matrix_set(model->V, model->K-1, 2, 0, 0.9711956316284261);
	matrix_set(model->V, model->K-1, 2, 1, 0.5010714686310176);

	// start test code //
	gensvm_write_model(model, "./data/test_write_model.txt");

	FILE *fid = fopen("./data/test_write_model.txt", "r");
	mu_assert(fid != NULL, "Couldn't open output file for reading");

	char buffer[GENSVM_MAX_LINE_LENGTH];

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "Output file for GenSVM (version 0.2.1)\n")
			== 0, "Line doesn't contain expected content (0).\n");

	// skip the time line
	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "\n") == 0,
		       	"Line doesn't contain expected content (1).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "Model:\n") == 0,
		       	"Line doesn't contain expected content (2).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "p = 0.9032800000000000\n") == 0,
		       	"Line doesn't contain expected content (3).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "lambda = 0.0130000000000000\n") == 0,
		       	"Line doesn't contain expected content (4).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "kappa = 1.1832000000000000\n") == 0,
		       	"Line doesn't contain expected content (5).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "epsilon = 1e-08\n") == 0,
		       	"Line doesn't contain expected content (6).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "weight_idx = 1\n") == 0,
		       	"Line doesn't contain expected content (7).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "\n") == 0,
		       	"Line doesn't contain expected content (8).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "Data:\n") == 0,
		       	"Line doesn't contain expected content (9).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "filename = ./data/test_file_read_data.txt\n")
			== 0, "Line doesn't contain expected content (10).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "n = 10\n") == 0,
		       	"Line doesn't contain expected content (11).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "m = 2\n") == 0,
		       	"Line doesn't contain expected content (12).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "K = 3\n") == 0,
		       	"Line doesn't contain expected content (13).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "\n") == 0,
		       	"Line doesn't contain expected content (14).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "Output:\n") == 0,
		       	"Line doesn't contain expected content (15).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "+0.4989893785603748 +0.0599082796573645\n")
		       	== 0, "Line doesn't contain expected content (16).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "+0.7918204761759593 +0.6456613497110559\n")
			== 0, "Line doesn't contain expected content (17).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "+0.9711956316284261 +0.5010714686310176\n")
			== 0, "Line doesn't contain expected content (18).\n");

	fclose(fid);

	// end test code //

	gensvm_free_model(model);

	return NULL;
}

char *test_gensvm_write_predictions()
{
	int n = 5,
	    m = 3,
	    K = 3;
	struct GenData *data = gensvm_init_data();
	long *predy = Calloc(long, n);

	data->n = n;
	data->m = m;
	data->K = K;

	data->Z = Calloc(double, data->n * (data->m+1));

	matrix_set(data->Z, data->m+1, 0, 0, 1.0);
	matrix_set(data->Z, data->m+1, 1, 0, 1.0);
	matrix_set(data->Z, data->m+1, 2, 0, 1.0);
	matrix_set(data->Z, data->m+1, 3, 0, 1.0);
	matrix_set(data->Z, data->m+1, 4, 0, 1.0);

	matrix_set(data->Z, data->m+1, 0, 1, 0.7065937536993949);
	matrix_set(data->Z, data->m+1, 0, 2, 0.7016517970438980);
	matrix_set(data->Z, data->m+1, 0, 3, 0.1548611397288129);
	matrix_set(data->Z, data->m+1, 1, 1, 0.4604987687863951);
	matrix_set(data->Z, data->m+1, 1, 2, 0.6374142980176117);
	matrix_set(data->Z, data->m+1, 1, 3, 0.0370930278245423);
	matrix_set(data->Z, data->m+1, 2, 1, 0.3798777132278375);
	matrix_set(data->Z, data->m+1, 2, 2, 0.5745070018747664);
	matrix_set(data->Z, data->m+1, 2, 3, 0.2570906697837264);
	matrix_set(data->Z, data->m+1, 3, 1, 0.2789376050039792);
	matrix_set(data->Z, data->m+1, 3, 2, 0.4853242744610165);
	matrix_set(data->Z, data->m+1, 3, 3, 0.1894010436762711);
	matrix_set(data->Z, data->m+1, 4, 1, 0.7630904372339489);
	matrix_set(data->Z, data->m+1, 4, 2, 0.1341546320318005);
	matrix_set(data->Z, data->m+1, 4, 3, 0.6827430912944857);

	predy[0] = 3;
	predy[1] = 2;
	predy[2] = 1;
	predy[3] = 2;
	predy[4] = 1;

	// start test code //

	gensvm_write_predictions(data, predy,
			"./data/test_write_predictions.txt");

	FILE *fid = fopen("./data/test_write_predictions.txt", "r");
	mu_assert(fid != NULL, "Couldn't open output file for reading");

	char buffer[GENSVM_MAX_LINE_LENGTH];
	// skip the first two lines
	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "5\n") == 0,
		       	"Line doesn't contain expected content (0).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "3\n") == 0,
		       	"Line doesn't contain expected content (1).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "0.7065937536993949 0.7016517970438980 "
				"0.1548611397288129 3\n") == 0,
		       	"Line doesn't contain expected content (2).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "0.4604987687863951 0.6374142980176117 "
				"0.0370930278245423 2\n") == 0,
		       	"Line doesn't contain expected content (3).\n");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "0.3798777132278375 0.5745070018747664 "
				"0.2570906697837264 1\n") == 0,
		       	"Line doesn't contain expected content (4).\n");


	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "0.2789376050039792 0.4853242744610165 "
				"0.1894010436762711 2\n") == 0,
		       	"Line doesn't contain expected content (5).\n");


	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "0.7630904372339489 0.1341546320318005 "
				"0.6827430912944857 1\n") == 0,
		       	"Line doesn't contain expected content (6).\n");

	fclose(fid);
	// end test code //

	gensvm_free_data(data);
	free(predy);

	return NULL;
}

char *test_gensvm_time_string()
{
	// not sure how to unit test this function.

	mu_test_missing();
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_gensvm_read_data);
	mu_run_test(test_gensvm_read_data_sparse);
	mu_run_test(test_gensvm_read_data_no_label);
	mu_run_test(test_gensvm_read_data_libsvm);
	mu_run_test(test_gensvm_read_data_libsvm_0based);
	mu_run_test(test_gensvm_read_data_libsvm_sparse);
	mu_run_test(test_gensvm_read_data_libsvm_no_label);

	mu_run_test(test_gensvm_read_model);
	mu_run_test(test_gensvm_write_model);
	mu_run_test(test_gensvm_write_predictions);
	mu_run_test(test_gensvm_time_string);

	return NULL;
}

RUN_TESTS(all_tests);
