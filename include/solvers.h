void solver_RK4( int FWTR, Vec M_diag, Mat K, Mat Ks, Vec f, 
	int nDOFreg, int *DOFx_reg, int *DOFy_reg, 
	double *node2xy, int *node2DOF, int node_load, 
	int nDOFsrf, int *DOFx_srf, int *DOFy_srf, 
	double h, double cp );