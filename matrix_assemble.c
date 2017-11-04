#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "petscmat.h"
#include "elements.h"

#include "assembly.h"

void element_reg( int elem, double lambda, double mu, double rho, double *node2xy, int *elem2node, double *k_elem, double *m_elem_diag, double *m_elem );

void assemble_matrix( int FWTR, Vec M_diag, Mat K, Mat Ks, int nElem, int *elem2loc, int nDOFsrf, int *node2DOF, int *elem2node, double *node2xy, double *cp ) {

	clock_t start = clock(), diff;
	int msec;

	int i0, i1, i2;

	double m_elem_diag[18], m_elem[324], k_elem[324];
	int address[18], address_srf[18];

	double lambda, mu, rho, cs;

	printf("----- matrix assembly -----\n");

	lambda = 1e8;
	mu = 1e8;
	rho = 2200.;
	cs = sqrt(lambda/rho);
	*cp = sqrt((lambda+2*mu)/rho);
	printf(" lambda=%e\n",lambda);
	printf("     mu=%e\n",mu);
	printf("    rho=%e\n",rho);
	printf("     cp=%9.4f\n",*cp);
	printf("     cs=%9.4f\n",cs);

	if( FWTR==0 ) {
		for (i0=0; i0<nElem; i0++) {

			element_reg( i0, lambda, mu, rho, node2xy, elem2node, k_elem, m_elem_diag, m_elem );

			for (i1=0; i1<2; i1++) {
				for (i2=0; i2<9; i2++) {
					address[i1*9+i2] = node2DOF[elem2node[i0*9+i2]*2+i1];
				}
			}

			VecSetValues( M_diag, 18, address, m_elem_diag, ADD_VALUES );
			MatSetValues( K, 18, address, 18, address, k_elem, ADD_VALUES );
		}

		VecAssemblyBegin( M_diag );
		VecAssemblyEnd  ( M_diag );
		MatAssemblyBegin( K, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd  ( K, MAT_FINAL_ASSEMBLY );

	} else if( FWTR==1 ) {
		for (i0=0; i0<nElem; i0++) {

			element_reg( i0, lambda, mu, rho, node2xy, elem2node, k_elem, m_elem_diag, m_elem );

			for (i1=0; i1<2; i1++) {
				for (i2=0; i2<9; i2++) {
					address    [i1*9+i2] = node2DOF[elem2node[i0*9+i2]*2+i1];
					address_srf[i1*9+i2] = -node2DOF[elem2node[i0*9+i2]*2+i1]-2;
				}
			}

			VecSetValues( M_diag, 18, address, m_elem_diag, ADD_VALUES );
			MatSetValues( K, 18, address, 18, address, k_elem, ADD_VALUES );
			if( elem2loc[i0]==0 ) { // for surface element
				// MatSetValues( Ms, 18, address, 18, address_srf, m_elem, ADD_VALUES );
				MatSetValues( Ks, 18, address, 18, address_srf, k_elem, ADD_VALUES );
			}
		}

		VecAssemblyBegin( M_diag );
		VecAssemblyEnd  ( M_diag );
		MatAssemblyBegin( K, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd  ( K, MAT_FINAL_ASSEMBLY );
		// MatAssemblyBegin( Ms, MAT_FINAL_ASSEMBLY );
		// MatAssemblyEnd  ( Ms, MAT_FINAL_ASSEMBLY );
		MatAssemblyBegin( Ks, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd  ( Ks, MAT_FINAL_ASSEMBLY );
	}

    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
	printf(" elapse=%4i seconds %4i milliseconds\n", msec/1000, msec%1000);
}