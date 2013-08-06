#include "libMSVMMaj.h"

/*
	Generate the simplex matrix. A pointer to the matrix
	must be given, and the matrix must have been
	allocated.
*/
void simplex_gen(long K, double *U)
{
	long i, j;
	for (i=0; i<K; i++) {
		for (j=0; j<K-1; j++) {
			if (i <= j) {
				matrix_set(U, K-1, i, j, -1.0/sqrt(2.0*(j+1)*(j+2)));
			} else if (i == j+1) {
				matrix_set(U, K-1, i, j, sqrt((j+1)/(2.0*(j+2))));
			} else {
				matrix_set(U, K-1, i, j, 0.0);
			}
		}
	}
}

/*
	Generate the category matrix R. The category matrix has 1's everywhere
	except at the column corresponding to the label of instance i.
*/
void category_matrix(struct Model *model, struct Data *dataset)
{
	long i, j;
	long n = model->n;
	long K = model->K;

	for (i=0; i<n; i++) {
		for (j=0; j<K; j++) {
			if (dataset->y[i] != j+1) {
				matrix_set(model->R, K, i, j, 1.0);
			}
		}
	}
}

void simplex_diff(struct Model *model, struct Data *data)
{
	long i, j, k;
	double value;

	long n = model->n;
	long K = model->K;

	for (i=0; i<n; i++) {
		for (j=0; j<K-1; j++) {
			for (k=0; k<K; k++) {
				value = matrix_get(model->U, K-1, data->y[i]-1, j);
				value -= matrix_get(model->U, K-1, k, j);
				matrix3_set(model->UU, K-1, K, i, j, k, value);
			}
		}
	}
}

/*
	Calculate the errors Q based on the current value of V.
	It is assumed that the memory for Q has already been allocated. 
	In addition, the matrix ZV is calculated here. It is assigned to a 
	pre-allocated block of memory, since it would be inefficient to keep
	reassigning this block at every iteration.
*/
void calculate_errors(struct Model *model, struct Data *data, double *ZV)
{
	long i, j, k;
	double a, value;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			n,
			K-1,
			m+1,
			1.0,
			data->Z,
			m+1,
			model->V,
			K-1,
			0.0,
			ZV,
			K-1);

	Memset(model->Q, double, n*K);
	for (i=0; i<n; i++) {
		for (j=0; j<K-1; j++) {
			a = matrix_get(ZV, K-1, i, j);
			for (k=0; k<K; k++) {
				value = a * matrix3_get(model->UU, K-1, K, i, j, k);
				matrix_add(model->Q, K, i, k, value);
			}
		}
	}
}	

/*
	Calculate the Huber hinge errors for each error in the matrix Q.
*/
void calculate_huber(struct Model *model) 
{
	long i, j;
	double q, value;

	for (i=0; i<model->n; i++) {
		for (j=0; j<model->K; j++) {
			q = matrix_get(model->Q, model->K, i, j);
			value = 0.0;
			if (q <= -model->kappa) {
				value = 1.0 - q - (model->kappa+1.0)/2.0;
			} else if (q <= 1.0) {
				value = 1.0/(2.0*model->kappa+2.0)*pow(1.0 - q, 2.0);
			}
			matrix_set(model->H, model->K, i, j, value);
		}
	}
}

/*
	Calculate the value of the loss function based on the current estimate
	of V. 
*/
double get_msvmmaj_loss(struct Model *model, struct Data *data, double *ZV)
{
	long i, j;
	long n = data->n;
	long K = data->K;
	long m = data->m;

	double value, rowvalue, loss = 0.0;
	
	calculate_errors(model, data, ZV);	
	calculate_huber(model);
	
	for (i=0; i<n; i++) {
		rowvalue = 0;
		value = 0;
		for (j=0; j<K; j++) {
			value = matrix_get(model->H, K, i, j);
			value = pow(value, model->p);
			value *= matrix_get(model->R, K, i, j);
			rowvalue += value;
		}
		rowvalue = pow(rowvalue, 1.0/(model->p));
		rowvalue *= model->rho[i];
		loss += rowvalue;
	}
	loss /= ((double) n);
	
	value = 0;
	for (i=1; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			value += pow(matrix_get(model->V, K-1, i, j), 2.0);
		}
	}
	loss += model->lambda * value;		

	return loss;
}

/*
	Training loop is defined here.
*/
void main_loop(struct Model *model, struct Data *data) 
{
	long i, j, it = 0;
	double L, Lbar;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	srand(time(NULL));

	double *B = Calloc(double, n*(K-1));
	double *ZV = Calloc(double, n*(K-1));
	double *ZAZ = Calloc(double, (m+1)*(m+1));	
	double *ZAZV = Calloc(double, (m+1)*(K-1));
	double *ZAZVT = Calloc(double, (m+1)*(K-1));

	info("Starting main loop.\n");
	info("Dataset:\n");
	info("\tn = %i\n", n);
	info("\tm = %i\n", m);
	info("\tK = %i\n", K);
	info("Parameters:\n");
	info("\tkappa = %f\n", model->kappa);
	info("\tp = %f\n", model->p);
	info("\tlambda = %15.16f\n", model->lambda);
	info("\tepsilon = %g\n", model->epsilon);
	info("\n");

	simplex_gen(model->K, model->U);
	simplex_diff(model, data);
	category_matrix(model, data);

	// Initialize V
	for (i=0; i<m+1; i++) 
		for (j=0; j<K-1; j++)
			matrix_set(model->V, K-1, i, j, 1.0);
			//matrix_set(model->V, K-1, i, j, -1.0+2.0*rnd());
	
	L = get_msvmmaj_loss(model, data, ZV);
	Lbar = L + 2.0*model->epsilon*L;

	while ((it < MAX_ITER) && (Lbar - L)/L > model->epsilon)
	{
		// ensure V contains newest V and Vbar contains V from previous
		msvmmaj_update(model, data, B, ZAZ,
				ZAZV, ZAZVT);
		if (it > 50)
			step_doubling(model);

		Lbar = L;
		L = get_msvmmaj_loss(model, data, ZV);

		if (it%100 == 0)
			info("iter = %li, L = %15.16f, Lbar = %15.16f, reldiff = %15.16f\n",
					it, L, Lbar, (Lbar - L)/L);
		it++;
	}

	info("optimization finished, iter = %li\n", it-1);

	for (i=0; i<K-1; i++)
		model->t[i] = matrix_get(model->V, K-1, 0, i);
	for (i=1; i<m+1; i++)
		for (j=0; j<K-1; j++)
			matrix_set(model->W, K-1, i-1, j, matrix_get(model->V, K-1, i, j));

	free(B);
	free(ZV);
	free(ZAZ);
	free(ZAZV);
	free(ZAZVT);
}

void step_doubling(struct Model *model)
{
	long i, j;

	long m = model->m;
	long K = model->K;

	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_mult(model->V, K-1, i, j, 2.0);
			matrix_add(model->V, K-1, i, j, -matrix_get(model->Vbar, K-1, i, j));
		}
	}
}

int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B, 
		int LDB)
{
	extern void dposv_(char *UPLO, int *Np, int *NRHSp, double *A, 
			int *LDAp, double *B, int *LDBp, int *INFOp);
	int INFO;
	dposv_(&UPLO, &N, &NRHS, A, &LDA, B, &LDB, &INFO);
	return INFO;
}

int  dsysv(char UPLO, int N, int NRHS, double *A, int LDA, int *IPIV,
		double *B, int LDB, double *WORK, int LWORK)
{
	extern void dsysv_(char *UPLO, int *Np, int *NRHSp, double *A, 
			int *LDAp, int *IPIV, double *B, int *LDBp, 
			double *WORK, int *LWORK, int *INFOp);
	int INFO;
	dsysv_(&UPLO, &N, &NRHS, A, &LDA, IPIV, B, &LDB, WORK, &LWORK, &INFO);
	return INFO;
}

void msvmmaj_update(struct Model *model, struct Data *data, 
		double *B, double *ZAZ, double *ZAZV, double *ZAZVT)
{
	// Because msvmmaj_update is always called after a call to 
	// get_msvmmaj_loss with the latest V, it is unnecessary to recalculate
	// the matrix ZV, the errors Q, or the Huber errors H. Awesome!
	int status, class;
	long i, j, k;
	double Avalue, Bvalue;
	double omega, value, a, b, q, h, r;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	double kappa = model->kappa;
	double p = model->p;
	double *rho = model->rho;

	const double a2g2 = 0.25*p*(2.0*p - 1.0)*pow((kappa+1.0)/2.0,p-2.0);
	const double in = 1.0/((double) n);

	Memset(B, double, n*(K-1));
	Memset(ZAZ, double, (m+1)*(m+1));
	b = 0;
	for (i=0; i<n; i++) {
		value = 0;
		omega = 0;
		for (j=0; j<K; j++) {
			h = matrix_get(model->H, K, i, j);
			r = matrix_get(model->R, K, i, j);
			value += (h*r > 0) ? 1 : 0;
			omega += pow(h, p)*r;
		}
		class = (value <= 1.0) ? 1 : 0;
		omega = (1.0/p)*pow(omega, 1.0/p - 1.0);

		Avalue = 0;
		if (class == 1) {
			for (j=0; j<K; j++) {
				q = matrix_get(model->Q, K, i, j);
				if (q <= -kappa) {
					a = 0.25/(0.5 - kappa/2.0 - q);
					b = 0.5;
				} else if (q <= 1.0) {
					a = 1.0/(2.0*kappa + 2.0);
					b = (1.0 - q)*a;
				} else {
					a = -0.25/(0.5 - kappa/2.0 - q);
					b = 0;
				}
				for (k=0; k<K-1; k++) {
					Bvalue = in*rho[i]*b*matrix3_get(model->UU, K-1, K, i, k, j);
					matrix_add(B, K-1, i, k, Bvalue);
				}
				Avalue += a*matrix_get(model->R, K, i, j);
			}
		} else {
			if (2.0 - p < 0.0001) {
				for (j=0; j<K; j++) {
					q = matrix_get(model->Q, K, i, j);
					if (q <= -kappa) {
						b = 0.5 - kappa/2.0 - q;
					} else if ( q <= 1.0) {
						b = pow(1.0 - q, 3.0)/(2.0*pow(kappa + 1.0, 2.0));
					} else {
						b = 0;
					}
					for (k=0; k<K-1; k++) {
						Bvalue = in*rho[i]*omega*b*matrix3_get(model->UU, K-1, K, i, k, j);
						matrix_add(B, K-1, i, k, Bvalue);
					}
				}
				Avalue = 1.5*(K - 1.0);
			} else {
				for (j=0; j<K; j++) {
					q = matrix_get(model->Q, K, i, j);
					if (q <= (p + kappa - 1.0)/(p - 2.0)) {
						a = 0.25*pow(p, 2.0)*pow(0.5 - kappa/2.0 - q, p - 2.0);
					} else if (q <= 1.0) {
						a = a2g2;
					} else {
						a = 0.25*pow(p, 2.0)*pow((p/(p - 2.0))*(0.5 - kappa/2.0 - q), p - 2.0);
						b = a*(2.0*q + kappa - 1.0)/(p - 2.0) + 0.5*p*pow((p/(p - 2.0))*(0.5 - kappa/2.0 - q), p - 1.0);
					}
					if (q <= -kappa) {
						b = 0.5*p*pow(0.5 - kappa/2.0 - q, p - 1.0);
					} else if ( q <= 1.0) {
						b = p*pow(1.0 - q, 2.0*p - 1.0)/pow(2*kappa+2.0, p);
					}
					for (k=0; k<K-1; k++) {
						Bvalue = in*rho[i]*omega*b*matrix3_get(model->UU, K-1, K, i, k, j);
						matrix_add(B, K-1, i, k, Bvalue);
					}
					Avalue += a*matrix_get(model->R, K, i, j);
				}
			}
			Avalue *= omega;
		}
		Avalue *= in * rho[i];

		// Now we calculate the matrix ZAZ. Since this is 
		// guaranteed to be symmetric, we only calculate the 
		// upper part of the matrix, and then copy this over
		// to the lower part after all calculations are done.
		// Note that the use of dsym is faster than dspr, even
		// though dspr uses less memory.
		cblas_dsyr(
				CblasRowMajor,
				CblasUpper,
				m+1,
				Avalue,
				&data->Z[i*(m+1)],
				1,
				ZAZ,
				m+1);
	}
	// Copy upper to lower (necessary because we need to switch
	// to Col-Major order for LAPACK).
	/*
	for (i=0; i<m+1; i++)
		for (j=0; j<m+1; j++)
			matrix_set(ZAZ, m+1, j, i, matrix_get(ZAZ, m+1, i, j));
	*/
	
	// Calculate the right hand side of the system we 
	// want to solve.
	cblas_dsymm(
			CblasRowMajor,
			CblasLeft,
			CblasUpper,
			m+1,
			K-1,
			1.0,
			ZAZ,
			m+1,
			model->V,
			K-1,
			0.0,
			ZAZV, 
			K-1);
	cblas_dgemm(
			CblasRowMajor,
			CblasTrans,
			CblasNoTrans,
			m+1,
			K-1,
			n,
			1.0,
			data->Z,
			m+1,
			B,
			K-1,
			1.0,
			ZAZV,
			K-1);
	/* 
	 * Add lambda to all diagonal elements except the 
	 * first one. 
	 */
	i = 0;
	for (j=0; j<m; j++)
		ZAZ[i+=m+1 + 1] += model->lambda;
	
	// For the LAPACK call we need to switch to Column-
	// Major order. This is unnecessary for the matrix
	// ZAZ because it is symmetric. The matrix ZAZV 
	// must be converted however.
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			ZAZVT[j*(m+1)+i] = ZAZV[i*(K-1)+j];
	
	// We use the lower ('L') part of the matrix ZAZ, 
	// because we have used the upper part in the BLAS
	// calls above in Row-major order, and Lapack uses
	// column major order. 
	status = dposv(
			'L',
			m+1,
			K-1,
			ZAZ,
			m+1,
			ZAZVT,
			m+1);

	if (status != 0) {
		// This step should not be necessary, as the matrix
		// ZAZ is positive semi-definite by definition. It 
		// is included for safety.
		fprintf(stderr, "Received nonzero status from dposv: %i\n", status);
		int *IPIV = malloc((m+1)*sizeof(int));
		double *WORK = malloc(1*sizeof(double));
		status = dsysv(
				'L',
				m+1,
				K-1,
				ZAZ,
				m+1,
				IPIV,
				ZAZVT,
				m+1,
				WORK,
				-1);
		WORK = (double *)realloc(WORK, WORK[0]*sizeof(double));
		status = dsysv(
				'L',
				m+1,
				K-1,
				ZAZ,
				m+1,
				IPIV,
				ZAZVT,
				m+1,
				WORK,
				sizeof(WORK)/sizeof(double));
		if (status != 0)
			fprintf(stderr, "Received nonzero status from dsysv: %i\n", status);
	}

	// Return to Row-major order. The matrix ZAZVT contains the 
	// solution after the dposv/dsysv call.
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			ZAZV[i*(K-1)+j] = ZAZVT[j*(m+1)+i];

	// Store the previous V in Vbar, assign the new V
	// (which is stored in ZAZVT) to the model, and give ZAZVT the
	// address of Vbar. This should ensure that we keep
	// re-using assigned memory instead of reallocating at every
	// update.
	/* See this answer: http://stackoverflow.com/q/13246615/
	 * For now we'll just do it by value until the rest is figured out.
	ptr = model->Vbar;
	model->Vbar = model->V;
	model->V = ZAZVT;
	ZAZVT = ptr;
	*/
	
	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_set(model->Vbar, K-1, i, j, matrix_get(model->V, K-1, i, j));
			matrix_set(model->V, K-1, i, j, matrix_get(ZAZV, K-1, i, j));
		}
	}
}

void initialize_weights(struct Data *data, struct Model *model)
{
	long *groups;
	long i;

	long n = model->n;
	long K = model->K;

	if (model->weight_idx == 1) {
		for (i=0; i<n; i++)
			model->rho[i] = 1.0;
	} 
	else if (model->weight_idx == 2) {
		groups = Calloc(long, K);
		for (i=0; i<n; i++) {
			groups[data->y[i]-1]++;
		}
		for (i=0; i<n; i++) {
			model->rho[i] = 1.0/((double) groups[data->y[i]-1]);
		}
	} else {
		fprintf(stderr, "Unknown weight specification.\n");
		exit(1);
	}

}

void predict_labels(struct Data *data, struct Model *model, long *predy)
{
	long i, j, k, label;
	double norm, min_dist;

	long n = data->n; // note that model->n is the size of the training sample.
	long m = data->m;
	long K = model->K; //data->K does not necessarily equal the original K.

	double *S = Calloc(double, K-1);
	double *ZV = Calloc(double, n*(K-1));
	double *U = Calloc(double, K*(K-1));

	// Get the simplex matrix
	simplex_gen(K, U);

	// Generate the simplex-space vectors
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			n,
			K-1,
			m+1,
			1.0,
			data->Z,
			m+1,
			model->V,
			K-1,
			0.0,
			ZV,
			K-1);

	// Calculate the distance to each of the vertices of the simplex.
	// The closest vertex defines the class label.
	for (i=0; i<n; i++) {
		label = 0;
		min_dist = 1000000000.0;
		for (j=0; j<K; j++) {
			for (k=0; k<K-1; k++) {
				S[k] = matrix_get(ZV, K-1, i, k) - matrix_get(U, K-1, j, k);
			}
			norm = cblas_dnrm2(K, S, 1);
			if (norm < min_dist) {
				label = j+1;
				min_dist = norm;
			}
		}
		predy[i] = label;
	}

	free(ZV);
	free(U);
	free(S);
}

double prediction_perf(struct Data *data, long *predy)
{
	long i, correct = 0;
	double performance;

	for (i=0; i<data->n; i++)
		if (data->y[i] == predy[i])
			correct++;

	performance = ((double) correct)/((double) data->n) * 100.0;

	return performance;
}
