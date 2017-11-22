#include <stdio.h>
#include <string.h>
#include "petscmat.h"
#include "petscviewer.h"

#include "model.h"
#include "assembly.h"
#include "solvers.h"

int main( int argc, char *args[] ) {

	int FWTR; // FW=0, TR=1
	int i0, i1;
	int nNode, nElem;
	double dim_x, dim_y;
	double *node2xy;
	int *elem2loc, *elem2node;
	int nDOF;
	int *node2DOF;
	int *DOFx, *DOFy;
	int nDOFsrf, *DOFx_srf, *DOFy_srf;
	int node_load;

	double h, cp, cs, dist_max;

	Vec M_diag;
	Mat Ks, K;
	Vec f;

PetscReal vecnorm;
PetscViewer viewer;

	PetscInitialize( &argc, &args, NULL, NULL );

	/* Is it FW or TR? */
	if( argc!=2 ) {
		printf(" input argument must be provided. (FW/TR)\n");
		PetscFinalize();
		return 0;
	} else if( !strcmp(args[1],"FW") ) {
		FWTR = 0;
		printf(" forward\n");
	} else if( !strcmp(args[1],"TR") ) {
		FWTR = 1;
		printf(" time reversal\n");
	} else {
		printf(" input argument must be provided. (FW/TR)\n");
		PetscFinalize();
		return 0;
	}

	/* Create model. */
	create_model( FWTR, &nNode, &dim_x, &dim_y, &node2xy, &nElem, &elem2node, &elem2loc, &nDOF, &node2DOF, &DOFx, &DOFy, &nDOFsrf, &DOFx_srf, &DOFy_srf, &h );

	/* Create PETSc objects. */
	VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, &M_diag );
	VecSetOption( M_diag, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );

	VecDuplicate( M_diag, &f );

	MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, nDOF, 72, NULL, 72, NULL, &K );
	MatSetOption( K, MAT_SYMMETRIC, PETSC_TRUE );
	MatSetOption( K, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE );

	MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, nDOFsrf, 72, NULL, 72, NULL, &Ks );
	MatSetOption( Ks, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE );

	/* Create matrices and force vector. */
	assemble_matrix( FWTR, M_diag, K, Ks, nElem, elem2loc, nDOFsrf, node2DOF, elem2node, node2xy, &cp, &cs );
	assemble_body_force( dim_x, dim_y, nNode, node2xy, node2DOF, f, &node_load, &dist_max );

// VecNorm( M_diag, NORM_2, &vecnorm );
// printf("norm(M_diag)=%f\n",vecnorm);
// MatNorm( K, NORM_1, &vecnorm );
// printf("norm(K     )=%f\n",vecnorm);
// MatNorm( Ms, NORM_1, &vecnorm );
// printf("norm(Ms    )=%f\n",vecnorm);
// MatNorm( Ks, NORM_1, &vecnorm );
// printf("norm(Ks    )=%f\n",vecnorm);
// VecNorm( f, NORM_2, &vecnorm );
// printf("norm(f     )=%f\n",vecnorm);
// exit(1);

// PetscViewerASCIIOpen( PETSC_COMM_WORLD, "Ms.m", &viewer );
// PetscViewerSetFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
// MatView( Ms, viewer );
// PetscViewerDestroy( &viewer );
// exit(1);
// // PetscViewerASCIIOpen( PETSC_COMM_WORLD, "M.m", &viewer );
// // PetscViewerSetFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
// // VecView( M_diag, viewer );
// // PetscViewerDestroy( &viewer );

	/* Solve. */
	solver_RK4( FWTR, M_diag, K, Ks, f, nDOF, DOFx, DOFy, node2xy, node2DOF, node_load, nDOFsrf, DOFx_srf, DOFy_srf, h, cp, cs, dist_max );

	/* Wrap up. */
	VecDestroy( &M_diag );
	MatDestroy( &K );
	// MatDestroy( &Ms );
	MatDestroy( &Ks );
	VecDestroy( &f );

// 	// printf("node\n");
// 	// for (i0 = 0; i0 < nNode; i0++) {
// 	// 	printf("%5.2f %5.2f", node2xy[2*i0], node2xy[2*i0 + 1]);
// 	// 	for (i1=0; i1<2; i1++) {
// 	// 		printf(" %3i ",node2DOF[i0*2+i1]);
// 	// 	}
// 	// 	printf("\n");
// 	// }

// 	// printf("element\n");
// 	// for( i0=0; i0<nElem; i0++){
// 	// 	for( i1=0; i1<9; i1++ ){
// 	// 		printf("%3i",elem2node[9*i0+i1]);
// 	// 	}
// 	// 	printf("\n");
// 	// }

	free( node2xy );
	free( elem2loc );
	free( elem2node );
	free( node2DOF );
	free( DOFx );
	free( DOFy );
	free( DOFx_srf );
	free( DOFy_srf );

	PetscFinalize();

	return 0;
}