
struct Task {
	enum KernelType kernel_type;
	int weight_idx;
	double epsilon;
	double p;
	double kappa;
	double lambda;
	double *kernel_param;
	struct MajData **data;
}


