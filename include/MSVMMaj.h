
#define MAX_ITER 10000000
#define MAX_LINE_LENGTH 1024

/*
	Model structure
*/
struct Model {
	int weight_idx;
	long K;
	long n;
	long m;
	double epsilon;
	double p;
	double kappa;
	double lambda;
	double *W;
	double *t;
	double *V;
	double *Vbar;
	double *U;
	double *UU;
	double *Q;
	double *H;
	double *R;
	double *rho;
	char *data_file;
};

/*
	Data structure
*/
struct Data {
	long K;
	long n;
	long m;
	long *y;
	double *Z;
};
