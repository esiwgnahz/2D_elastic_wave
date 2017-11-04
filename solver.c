#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "petscmat.h"
#include "petscviewerhdf5.h"
#include "petscis.h"

double loading_fw_time_signal( double t, double offset );

double dabsmax( int nsize, const double *arr ) {
	
	int i0;
	double maximum;

	maximum = 0.;

	for( i0=0; i0<nsize; i0++ ) {
		if( maximum<fabs(arr[i0]) ) maximum = arr[i0];
	}

	return maximum;
}

void solver_RK4( int FWTR, Vec M_diag, Mat K, Mat Ks, Vec f, int nDOFreg, int *DOFx_reg, int *DOFy_reg, double *node2xy, int *node2DOF, int node_load, int nDOFsrf, int *DOFx_srf, int *DOFy_srf ) {

	clock_t start = clock(), diff;
	int msec;

	Vec k11, k12, k21, k22, k31, k32, k41, k42;
	Vec x1n, x2n;
	Vec temp;

	double maxval_x, maxval_x_new;
	double maxval_y, maxval_y_new;

	int i0, i1;
	int nTstep;
	double dt;
	double offset;

	// print load/target response
	PetscInt DOF_load[2];
	double disp_load[2];

	// print entire domain response
	IS isx_reg, isy_reg;
	Vec ux_reg, uy_reg;
	const double *ux_reg_ftr, *uy_reg_ftr;
	int timestep;

	// print surface response
	IS isx_srf, isy_srf;
	Vec ux_srf, uy_srf, u_srf;
	double *ux_srf_ftr, *uy_srf_ftr;

	// read final state of forward step
	double *ux_fin_ftr, *uy_fin_ftr;

	// files
	FILE *fid_ux_srf, *fid_uy_srf;
	FILE *fid_ux_reg, *fid_uy_reg;
	FILE *fid_ux_trg, *fid_uy_trg;
	FILE *fid_ux_fin, *fid_uy_fin;

	printf("----- solver RK4 -----\n");

	nTstep = 425;
	dt = .001;
	offset = .01;

	DOF_load[0] = node2DOF[node_load*2  ];
	DOF_load[1] = node2DOF[node_load*2+1];

	VecDuplicate( M_diag, &k11 ); VecDuplicate( M_diag, &k12 );
	VecDuplicate( M_diag, &k21 ); VecDuplicate( M_diag, &k22 );
	VecDuplicate( M_diag, &k31 ); VecDuplicate( M_diag, &k32 );
	VecDuplicate( M_diag, &k41 ); VecDuplicate( M_diag, &k42 );
	VecDuplicate( M_diag, &x1n ); VecDuplicate( M_diag, &x2n );
	VecDuplicate( M_diag, &temp );

	/* Vec to read surface response of forward step. */
	if( FWTR==1 ) {
		VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOFsrf, &u_srf );
		ux_srf_ftr = (double*)malloc(nDOFsrf/2*sizeof(double));
		uy_srf_ftr = (double*)malloc(nDOFsrf/2*sizeof(double));
	}

	VecReciprocal( M_diag );
	MatDiagonalScale( K, M_diag, NULL );
	MatScale( K, -1. );
	if( FWTR==0 ){
		VecPointwiseMult( f, f, M_diag );
	} else if( FWTR==1 ) {
		MatScale( Ks, -1. );
	}

	VecSet( x1n, 0. );
	VecSet( x2n, 0. );
	maxval_x = 0.; maxval_y = 0.;
	timestep = 0;
	if( FWTR==1 ) VecSet( u_srf, 0. );

	ISCreateGeneral( PETSC_COMM_WORLD, nDOFreg/2, DOFx_reg, PETSC_COPY_VALUES, &isx_reg );
	ISCreateGeneral( PETSC_COMM_WORLD, nDOFreg/2, DOFy_reg, PETSC_COPY_VALUES, &isy_reg );

	ISCreateGeneral( PETSC_COMM_WORLD, nDOFsrf/2, DOFx_srf, PETSC_COPY_VALUES, &isx_srf );
	ISCreateGeneral( PETSC_COMM_WORLD, nDOFsrf/2, DOFy_srf, PETSC_COPY_VALUES, &isy_srf );

	// /* Read final state of forward step. */
	// if( FWTR==1 ) {
	// 	ux_fin_ftr = (double*)malloc(nDOFreg/2*sizeof(double));
	// 	uy_fin_ftr = (double*)malloc(nDOFreg/2*sizeof(double));

	// 	fid_ux_fin = fopen("output/ux_reg.dat","rb");
	// 	fid_uy_fin = fopen("output/uy_reg.dat","rb");

	// 	fseek( fid_ux_fin, -(nDOFreg/2)*sizeof(double), SEEK_END );
	// 	fseek( fid_uy_fin, -(nDOFreg/2)*sizeof(double), SEEK_END );

	// 	fread( ux_fin_ftr, sizeof(double), nDOFreg/2, fid_ux_fin );
	// 	fread( uy_fin_ftr, sizeof(double), nDOFreg/2, fid_uy_fin );

	// 	VecSetValues( x1n, nDOFreg/2, DOFx_reg, ux_fin_ftr, INSERT_VALUES );
	// 	VecSetValues( x1n, nDOFreg/2, DOFy_reg, uy_fin_ftr, INSERT_VALUES );

	// 	VecAssemblyBegin( x1n );
	// 	VecAssemblyEnd  ( x1n );

	// 	fclose( fid_ux_fin );
	// 	fclose( fid_uy_fin );

	// 	free( ux_fin_ftr );
	// 	free( uy_fin_ftr );
	// }

	/* Open file for input and output. */
	fid_ux_reg = fopen("output/ux_reg.dat","wb");
	fid_uy_reg = fopen("output/uy_reg.dat","wb");
	if( FWTR==0 ) {
		fid_ux_srf = fopen("output/ux_srf.dat","wb");
		fid_uy_srf = fopen("output/uy_srf.dat","wb");
	} else if( FWTR==1 ) {
		fid_ux_srf = fopen("output/ux_srf.dat","rb");
		fid_uy_srf = fopen("output/uy_srf.dat","rb");
	}
	fid_ux_trg = fopen("output/ux_trg.dat","wb");
	fid_uy_trg = fopen("output/uy_trg.dat","wb");

	/* Start time stepping. */
	for( i0=0; i0<nTstep; i0++ ) {

		/* loading */
		VecSet( k12, 0. );
		VecSet( k22, 0. );
		VecSet( k32, 0. );
		VecSet( k42, 0. );
		if( FWTR==0 ) {
			VecAXPY( k12, loading_fw_time_signal( (       i0   )*dt, offset ), f );
			VecAXPY( k22, loading_fw_time_signal( (       i0+.5)*dt, offset ), f );
			VecAXPY( k32, loading_fw_time_signal( (       i0+.5)*dt, offset ), f );
			VecAXPY( k42, loading_fw_time_signal( (       i0+1.)*dt, offset ), f );
		} else if( FWTR==1 ) {
			// VecAXPY( k12, loading_fw_time_signal( (nTstep-i0   )*dt, offset ), f );
			// VecAXPY( k22, loading_fw_time_signal( (nTstep-i0-.5)*dt, offset ), f );
			// VecAXPY( k32, loading_fw_time_signal( (nTstep-i0-.5)*dt, offset ), f );
			// VecAXPY( k42, loading_fw_time_signal( (nTstep-i0-1.)*dt, offset ), f );

			fseek( fid_ux_srf, -(i0+1)*(nDOFsrf/2)*sizeof(double), SEEK_END );
			fread( ux_srf_ftr, sizeof(double), nDOFsrf/2, fid_ux_srf );
			VecSetValues( u_srf, nDOFsrf/2, DOFx_srf, ux_srf_ftr, INSERT_VALUES );

			fseek( fid_uy_srf, -(i0+1)*(nDOFsrf/2)*sizeof(double), SEEK_END );
			fread( uy_srf_ftr, sizeof(double), nDOFsrf/2, fid_uy_srf );
			VecSetValues( u_srf, nDOFsrf/2, DOFy_srf, uy_srf_ftr, INSERT_VALUES );

			VecAssemblyBegin( u_srf );
			VecAssemblyEnd  ( u_srf );

			MatMultAdd( Ks, u_srf, k12, k12 );
			MatMultAdd( Ks, u_srf, k22, k22 );
			MatMultAdd( Ks, u_srf, k32, k32 );
			MatMultAdd( Ks, u_srf, k42, k42 );
			VecPointwiseMult( k12, M_diag, k12 );
			VecPointwiseMult( k22, M_diag, k22 );
			VecPointwiseMult( k32, M_diag, k32 );
			VecPointwiseMult( k42, M_diag, k42 );
		}

		/* Runge-Kutta 4 */
		VecCopy( x2n, k11 );
		MatMultAdd( K, x1n, k12, k12 );

		VecWAXPY( k21, .5*dt, k12, x2n );
		VecWAXPY( temp, .5*dt, k11, x1n );
		MatMultAdd( K, temp, k22, k22 );

		VecWAXPY( k31, .5*dt, k22, x2n );
		VecWAXPY( temp, .5*dt, k21, x1n );
		MatMultAdd( K, temp, k32, k32 );

		VecWAXPY( k41, dt, k32, x2n );
		VecWAXPY( temp, dt, k31, x1n );
		MatMultAdd( K, temp, k42, k42 );

		/* Update solution vectors. */
		VecSet( temp, 0. );
		VecAXPY( temp, 1., k11 );
		VecAXPY( temp, 2., k21 );
		VecAXPY( temp, 2., k31 );
		VecAXPY( temp, 1., k41 );
		VecAXPY( x1n, dt/6., temp );

		VecSet( temp, 0. );
		VecAXPY( temp, 1., k12 );
		VecAXPY( temp, 2., k22 );
		VecAXPY( temp, 2., k32 );
		VecAXPY( temp, 1., k42 );
		VecAXPY( x2n, dt/6., temp );

		/* Print surface response. */
		if( FWTR==0 ) {
			VecGetSubVector( x1n, isx_srf, &ux_srf );
			VecGetSubVector( x1n, isy_srf, &uy_srf );

			VecGetArrayRead( ux_srf, &ux_srf_ftr );
			VecGetArrayRead( uy_srf, &uy_srf_ftr );

			fwrite( ux_srf_ftr, sizeof(PetscScalar), nDOFsrf/2, fid_ux_srf );
			fwrite( uy_srf_ftr, sizeof(PetscScalar), nDOFsrf/2, fid_uy_srf );

			VecRestoreArrayRead( ux_srf, &ux_srf_ftr );
			VecRestoreArrayRead( uy_srf, &uy_srf_ftr );

			VecRestoreSubVector( x1n, isx_srf, &ux_srf );
			VecRestoreSubVector( x1n, isy_srf, &uy_srf );
		}

		/* Print entire response. */
		if( i0%5==0 ) {
			VecGetSubVector( x1n, isx_reg, &ux_reg );
			VecGetSubVector( x1n, isy_reg, &uy_reg );

			VecGetArrayRead( ux_reg, &ux_reg_ftr );
			VecGetArrayRead( uy_reg, &uy_reg_ftr );

			fwrite( ux_reg_ftr, sizeof(PetscScalar), nDOFreg/2, fid_ux_reg );
			fwrite( uy_reg_ftr, sizeof(PetscScalar), nDOFreg/2, fid_uy_reg );

			maxval_x_new = dabsmax( nDOFreg/2, ux_reg_ftr );
			if( maxval_x<maxval_x_new ) maxval_x = maxval_x_new;
			maxval_y_new = dabsmax( nDOFreg/2, uy_reg_ftr );
			if( maxval_y<maxval_y_new ) maxval_y = maxval_y_new;

			VecRestoreArrayRead( ux_reg, &ux_reg_ftr );
			VecRestoreArrayRead( uy_reg, &uy_reg_ftr );

			VecRestoreSubVector( x1n, isx_reg, &ux_reg );
			VecRestoreSubVector( x1n, isy_reg, &uy_reg );

			// timing
			diff = clock() - start;
    		msec = diff * 1000 / CLOCKS_PER_SEC;
			// printf(" step %i/%i elapse=%4i seconds %4i milliseconds\n", i0+1,nTstep, msec/1000, msec%1000);
			start = clock();
		}

		/* Print response on the loading nodes. */
		VecGetValues( x1n, 2, DOF_load, disp_load );
		fwrite( &disp_load[0], sizeof(double), 1, fid_ux_trg );
		fwrite( &disp_load[1], sizeof(double), 1, fid_uy_trg );
		printf("t=%5.3f, x=%5.2f, y=%5.2f, ux_trg=%11.4e, uy_trg=%11.4e\n",(i0+1)*dt,node2xy[node_load*2],node2xy[node_load*2+1],disp_load[0],disp_load[1]);
	}

	ISDestroy( &isx_reg );
	ISDestroy( &isy_reg );
	ISDestroy( &isx_srf );
	ISDestroy( &isy_srf );

	if( FWTR==1 ) {
		free( ux_srf_ftr );
		free( uy_srf_ftr );
	}

	fclose( fid_ux_reg );
	fclose( fid_uy_reg );
	fclose( fid_ux_srf );
	fclose( fid_uy_srf );
	fclose( fid_ux_trg );
	fclose( fid_uy_trg );

	printf("maxval x=%e\n",maxval_x);
	printf("maxval y=%e\n",maxval_y);
}